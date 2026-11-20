from pathlib import Path
import anndata as ad
import pandas as pd
import torch

from config import SCGPT_DIR, GENEFORMER_DIR, GENEFORMER_ASSETS, OUTPUT_DIR
from data import SingleCellDataset
from models import MODELS


def load_model(name, device):
    if name == "scgpt":
        return MODELS[name](SCGPT_DIR, device)
    return MODELS[name](GENEFORMER_DIR, GENEFORMER_ASSETS, device)


def run(model_name, datasets, batch_size=None, split="all", output_dir=OUTPUT_DIR):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = load_model(model_name, device)  # Once per model, not once per dataset.
    batch_size = batch_size or (32 if model_name == "scgpt" else 8)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    for name in datasets:
        data = SingleCellDataset(name, split=split)
        try:
            embeddings = model.embed(data, batch_size=batch_size)
            result = ad.AnnData(obs=pd.DataFrame({"label": data.labels}, index=data.cell_ids))
            for layer in range(embeddings.shape[1]):
                result.obsm[f"X_{model_name}_layer_{layer + 1}"] = embeddings[:, layer, :]
            target = Path(output_dir) / f"{name}_{split}_{model_name}.h5ad"
            result.write_h5ad(target)
            print(f"Saved {target}: {embeddings.shape}")
        finally:
            data.close()
