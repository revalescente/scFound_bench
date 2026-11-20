import pickle
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from torch.nn.utils.rnn import pad_sequence
from tqdm.auto import tqdm
from transformers import BertForMaskedLM


class GeneformerCells(Dataset):
    def __init__(self, data, tokens, medians, mapping, max_length):
        self.data, self.tokens, self.max_length = data, tokens, max_length
        genes = data.ensembl_ids
        if genes is None:
            genes = data.gene_names
        
        genes = [mapping.get(g, g) for g in genes]
        self.positions = np.array([i for i, g in enumerate(genes) if g in tokens and g in medians and medians[g] > 0])
        selected = [genes[i] for i in self.positions]
        unique, self.inverse = np.unique(selected, return_inverse=True)
        self.gene_tokens = np.array([tokens[g] for g in unique])
        self.gene_medians = np.array([medians[g] for g in unique])
        print(f"Geneformer genes matched: {len(selected)}/{len(genes)}")

    def __len__(self):
        return len(self.data)

    def __getitem__(self, i):
        counts = self.data[i][self.positions]
        merged = np.bincount(self.inverse, weights=counts, minlength=len(self.gene_tokens))
        total_counts = merged.sum()
        normalized = (
            merged / total_counts * 10000
        ) / self.gene_medians
        
        scores = normalized
        expressed = np.flatnonzero(scores > 0)
        order = expressed[np.argsort(-scores[expressed], kind="stable")]
        genes = self.gene_tokens[order[:self.max_length - 2]]
        return [self.tokens["<cls>"], *genes.tolist(), self.tokens["<eos>"]]


def make_collator(pad_id):
    def collate(batch):
        lengths = torch.tensor([len(row) for row in batch])
        ids = pad_sequence([torch.tensor(row) for row in batch], batch_first=True, padding_value=pad_id)
        return ids, lengths
    return collate


class Geneformer:
    def __init__(self, path, assets, device, max_length=4096):
        self.device, self.max_length = torch.device(device), max_length
        names = ("token_dictionary_gc104M.pkl", "gene_median_dictionary_gc104M.pkl", "ensembl_mapping_dict_gc104M.pkl")
        self.tokens, self.medians, self.mapping = [pickle.loads((assets / name).read_bytes()) for name in names]
        self.model = BertForMaskedLM.from_pretrained(
            path, local_files_only=True, output_hidden_states=True
        ).to(self.device).eval()

    def embed(self, data, batch_size=8):
        cells = GeneformerCells(data, self.tokens, self.medians, self.mapping, self.max_length)
        loader = DataLoader(cells, batch_size=batch_size, num_workers=4,
                            pin_memory=self.device.type == "cuda",
                            collate_fn=make_collator(self.tokens["<pad>"]))
        collected = []
        with torch.inference_mode():
            for ids, lengths in tqdm(loader, desc="Geneformer"):
                ids = ids.to(self.device, non_blocking=True)
                lengths = lengths.to(self.device)
                positions = torch.arange(ids.size(1), device=self.device)[None, :]
                valid = (positions > 0) & (positions < lengths[:, None] - 1)
                mask = (positions < lengths[:, None]).long()
                with torch.autocast(device_type=self.device.type, enabled=self.device.type == "cuda"):
                    hidden = self.model.bert(input_ids=ids, attention_mask=mask, output_hidden_states=True).hidden_states[1:]
                    embeddings = torch.stack([
                        (h * valid.unsqueeze(-1)).sum(1) / valid.sum(1).clamp(min=1).unsqueeze(-1)
                        for h in hidden
                    ], dim=1)
                collected.append(embeddings.float().cpu())
        return torch.cat(collected).numpy()
