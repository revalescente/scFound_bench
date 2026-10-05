#!/usr/bin/env python3
"""scGPT: estrazione embeddings per layer, un CSV per layer per dataset."""
import json
from pathlib import Path
import torch
from scgpt.model import TransformerModel
from scgpt.tokenizer import GeneVocab
from scgpt.data_collator import DataCollator
from torch.utils.data import DataLoader, SequentialSampler
import anndata as an
import numpy as np
import os
from tqdm import tqdm
import pandas as pd

# ── Costanti ───────────────────────────────────────────────────────────────
MODEL_DIR         = Path("/projects/shared/intronic_bam/scGPT_supplements/save/scGPT_human")
OUTPUT_BASE_DIR   = Path("/projects/shared/intronic_bam/datasets/embeddings/scgpt_layers")
MODEL_CONFIG_FILE = MODEL_DIR / "args.json"
MODEL_FILE        = MODEL_DIR / "best_model.pt"
VOCAB_FILE        = MODEL_DIR / "vocab.json"   # era model_dir (bug)

BATCH_SIZE = 8    # era undefined
MAX_LENGTH = 1200

DATASETS_CONFIG = [
    {
        "name": "be1",
        "h5ad_path": "/projects/shared/intronic_bam/datasets/anndata/be1.h5ad",
        "meta_cols": ["Barcode", "Sample"],
        "gene_col": "gene_names",
        "layer": "counts",   # mancava la virgola prima (bug)
    }
]

# ── Vocab ──────────────────────────────────────────────────────────────────
vocab = GeneVocab.from_file(VOCAB_FILE)
for s in ["<pad>", "<cls>", "<eoc>", "<mask>", "<unk>"]:
    if s not in vocab:
        vocab.append_token(s)
vocab.set_default_index(vocab["<pad>"])

# ── Configurazione e modello ───────────────────────────────────────────────
with open(MODEL_CONFIG_FILE) as f:
    model_configs = json.load(f)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model = TransformerModel(
    ntoken=len(vocab),
    d_model=model_configs["embsize"],
    nhead=model_configs["nheads"],
    d_hid=model_configs["d_hid"],
    nlayers=model_configs["nlayers"],
    nlayers_cls=model_configs.get("n_layers_cls", 3),
    vocab=vocab,
    pad_token="<pad>",
    pad_value=model_configs.get("pad_value", 0),
    do_mvc=model_configs.get("GEPC", False),
)

model_dict = model.state_dict()
pretrained_dict = torch.load(MODEL_FILE, map_location=device)
pretrained_dict = {
    k: v for k, v in pretrained_dict.items()
    if k in model_dict and v.shape == model_dict[k].shape
}
model_dict.update(pretrained_dict)
model.load_state_dict(model_dict)
model.to(device)
model.eval()
print(f"Modello caricato su: {device}")


# ── Dataset ────────────────────────────────────────────────────────────────
class GeneDataset(torch.utils.data.Dataset):
    def __init__(self, count_matrix, gene_ids):
        self.count_matrix = count_matrix
        self.gene_ids = gene_ids

    def __len__(self):
        return len(self.count_matrix)

    def __getitem__(self, idx):
        row = self.count_matrix[idx]
        nz = np.nonzero(row)[0]
        genes  = np.insert(self.gene_ids[nz], 0, vocab["<cls>"])
        values = np.insert(row[nz], 0, model_configs["pad_value"])
        return {
            "id": idx,
            "genes": torch.from_numpy(genes).long(),
            "expressions": torch.from_numpy(values).float(),
        }


# ── Estrazione embeddings ──────────────────────────────────────────────────
def get_all_layer_embeddings(model, data_loader, vocab, model_configs, normalize=True):
    device = next(model.parameters()).device
    emb_layers, hooks = {}, []

    def get_layer_embed(name):
        def hook(module, inp, out):
            out = out[0] if isinstance(out, (tuple, list)) else out
            if isinstance(out, torch.Tensor):
                # estrai subito solo il CLS token → [bs, 512] invece di [bs, 1200, 512]
                emb_layers[name] = out[:, 0, :].detach().cpu()
        return hook

    for name, module in model.named_modules():
        if "transformer_encoder.layers." in name and name.split(".")[-1].isdigit():
            hooks.append(module.register_forward_hook(get_layer_embed(name)))

    collected = {}
    try:
        with torch.no_grad(), torch.amp.autocast("cuda", enabled=True):
            for data_dict in tqdm(data_loader, desc="Embedding cells"):
                emb_layers.clear()
                input_gene_ids = data_dict["gene"].to(device)
                src_key_padding_mask = input_gene_ids.eq(vocab[model_configs["pad_token"]])
                model._encode(
                    input_gene_ids,
                    data_dict["expr"].to(device),
                    src_key_padding_mask=src_key_padding_mask,
                )
                for name, tensor in emb_layers.items():
                    collected.setdefault(name, []).append(tensor.numpy())  # già [bs, 512]
    finally:
        for h in hooks:
            h.remove()

    layer_names = sorted(collected, key=lambda n: int(n.split(".")[-1]))
    all_embs = np.stack([np.concatenate(collected[n], axis=0) for n in layer_names], axis=1)

    if normalize:
        norms = np.linalg.norm(all_embs, axis=2, keepdims=True)
        all_embs /= np.maximum(norms, 1e-8)

    return all_embs, layer_names


# ── Loop sui dataset ───────────────────────────────────────────────────────
for cfg in DATASETS_CONFIG:
    name = cfg["name"]
    print(f"\n=== {name} ===")

    save_dir = OUTPUT_BASE_DIR / name
    save_dir.mkdir(parents=True, exist_ok=True)

    adata = an.read_h5ad(cfg["h5ad_path"])
    adata.X = adata.layers[cfg["layer"]].copy()

    gene_col = cfg["gene_col"]
    adata.var["id_in_vocab"] = [vocab[g] if g in vocab else -1 for g in adata.var[gene_col]]
    adata = adata[:, adata.var["id_in_vocab"] >= 0].copy()  # .copy() evita view issues

    gene_ids = np.array(vocab(adata.var[gene_col].tolist()), dtype=int)
    count_matrix = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()

    dataset = GeneDataset(count_matrix, gene_ids)
    collator = DataCollator(
        do_padding=True,
        pad_token_id=vocab[model_configs["pad_token"]],
        pad_value=model_configs["pad_value"],
        do_mlm=False,
        do_binning=True,
        max_length=MAX_LENGTH,
        sampling=True,
        keep_first_n_tokens=1,
    )
    data_loader = DataLoader(
        dataset,
        batch_size=BATCH_SIZE,
        sampler=SequentialSampler(dataset),
        collate_fn=collator,
        drop_last=False,
        num_workers=min(len(os.sched_getaffinity(0)), BATCH_SIZE),
        pin_memory=True,
    )

    embs, layer_names = get_all_layer_embeddings(model, data_loader, vocab, model_configs)
    print(f"Shape embeddings: {embs.shape}")  # (n_cells, 12, 512)

    for i, lname in enumerate(layer_names):
        layer_idx = lname.split(".")[-1]
        df = pd.DataFrame(embs[:, i, :], index=adata.obs.index)
        df.index.name = "Barcode"
        df.insert(0, "Sample", adata.obs["Sample"].values)
        df.to_csv(save_dir / f"be1_scgpt_layer_{layer_idx}.csv")

    print(f"Salvati {len(layer_names)} CSV in {save_dir}")
