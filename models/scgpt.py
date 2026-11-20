import json
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm.auto import tqdm

from scgpt.data_collator import DataCollator
from scgpt.model import TransformerModel
from scgpt.tokenizer import GeneVocab
from scgpt.utils import load_pretrained


class ScGPTCells(Dataset):
    def __init__(self, data, vocab, pad_value):
        self.data, self.vocab, self.pad_value = data, vocab, pad_value
        ids = np.array([vocab[g] if g in vocab else -1 for g in data.gene_names])
        self.valid = ids >= 0
        self.gene_ids = ids[self.valid]
        print(f"scGPT genes matched: {self.valid.sum()}/{len(ids)}")

    def __len__(self):
        return len(self.data)

    def __getitem__(self, i):
        values = self.data[i][self.valid]
        used = np.flatnonzero(values)
        genes = np.r_[self.vocab["<cls>"], self.gene_ids[used]]
        expr = np.r_[self.pad_value, values[used]]
        return {"id": i, "genes": torch.as_tensor(genes, dtype=torch.long),
                "expressions": torch.as_tensor(expr, dtype=torch.float32)}


class ScGPT:
    def __init__(self, path, device, max_length=1200):
        self.device, self.max_length = torch.device(device), max_length
        vocab = GeneVocab.from_file(path / "vocab.json")
        for token in ("<pad>", "<cls>", "<eoc>"):
            if token not in vocab:
                vocab.append_token(token)
        vocab.set_default_index(vocab["<pad>"])
        self.vocab = vocab
        cfg = json.loads((path / "args.json").read_text())
        self.cfg = cfg
        self.model = TransformerModel(
            ntoken=len(vocab), d_model=cfg["embsize"], nhead=cfg["nheads"],
            d_hid=cfg["d_hid"], nlayers=cfg["nlayers"],
            nlayers_cls=cfg["n_layers_cls"], n_cls=1, vocab=vocab,
            dropout=cfg["dropout"], pad_token=cfg["pad_token"],
            pad_value=cfg["pad_value"], do_mvc=True, do_dab=False,
            use_batch_labels=False, domain_spec_batchnorm=False,
            explicit_zero_prob=False, use_fast_transformer=False,
            fast_transformer_backend="flash", pre_norm=False,
        )
        load_pretrained(self.model, torch.load(path / "best_model.pt", map_location=self.device), verbose=False)
        self.model.to(self.device).eval()

    def embed(self, data, batch_size=32):
        cells = ScGPTCells(data, self.vocab, self.cfg["pad_value"])
        collate = DataCollator(
            do_padding=True, pad_token_id=self.vocab[self.cfg["pad_token"]],
            pad_value=self.cfg["pad_value"], do_mlm=False, do_binning=True,
            max_length=self.max_length, sampling=True, keep_first_n_tokens=1,
        )
        loader = DataLoader(cells, batch_size=batch_size, collate_fn=collate, num_workers=0,
                            pin_memory=self.device.type == "cuda")
        collected, layer_outputs, hooks = [], [], []

        def hook(module, args, output):
            layer_outputs.append(output[:, 0, :])

        for layer in self.model.transformer_encoder.layers:
            hooks.append(layer.register_forward_hook(hook))
        try:
            with torch.inference_mode():
                for batch in tqdm(loader, desc="scGPT"):
                    genes = batch["gene"].to(self.device, non_blocking=True)
                    expr = batch["expr"].to(self.device, non_blocking=True)
                    layer_outputs.clear()
                    with torch.autocast(device_type=self.device.type, enabled=self.device.type == "cuda"):
                        self.model._encode(
                            genes, expr,
                            src_key_padding_mask=genes.eq(self.vocab[self.cfg["pad_token"]]),
                            batch_labels=None,
                        )
                    collected.append(torch.stack(layer_outputs, dim=1).float().cpu())
        finally:
            for h in hooks:
                h.remove()

        embeddings = torch.cat(collected).numpy()
        return embeddings / np.maximum(np.linalg.norm(embeddings, axis=-1, keepdims=True), 1e-12)
