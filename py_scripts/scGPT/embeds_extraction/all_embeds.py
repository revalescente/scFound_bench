#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scgpt as scg


MODEL_DIR = Path("/projects/shared/intronic_bam/scGPT_supplements/save/scGPT_human")
OUTPUT_BASE_DIR = Path("/projects/shared/intronic_bam/datasets/embeddings/scgpt")

DATASETS_CONFIG = [
    {
        "name": "be1",
        "h5ad_path": "/projects/shared/intronic_bam/datasets/anndata/be1.h5ad",
        "meta_cols": ["Barcode", "Sample"],
        "gene_col": "gene_names",
        "preprocess": "normalize",  # needs normalize_total + log1p
        "layer": "counts",
    },
    {
        "name": "cb",
        "h5ad_path": "/projects/shared/intronic_bam/datasets/anndata/cord_blood.h5ad",
        "meta_cols": ["celltype"],
        "gene_col": "gene_names",
        "preprocess": "normalize",
        "layer": "counts",
    },
    {
        "name": "mix",
        "h5ad_path": "/projects/shared/intronic_bam/datasets/anndata/sc_mix.h5ad",
        "meta_cols": ["cell_line"],
        "gene_col": "gene_names",
        "preprocess": "logcounts",  # already has logcounts layer
        "layer": "logcounts",
    },
]


def preprocess(adata: sc.AnnData, config: dict) -> sc.AnnData:
    if config["preprocess"] == "normalize":
        adata.X = adata.layers[config["layer"]].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif config["preprocess"] == "logcounts":
        adata.X = adata.layers[config["layer"]].copy()
        if "gene_names" not in adata.var.columns:
            adata.var["gene_names"] = adata.var_names
    return adata


def extract_embeddings(config: dict, batch_size: int = 64) -> None:
    dataset_name = config["name"]

    print(f"\n==========================================")
    print(f" Extracting Embeddings for: {dataset_name.upper()}")
    print(f"==========================================")

    output_dir = OUTPUT_BASE_DIR / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(config["h5ad_path"])
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    adata = preprocess(adata, config)

    embedded = scg.tasks.embed_data(
        adata,
        MODEL_DIR,
        gene_col=config["gene_col"],
        batch_size=batch_size,
    )

    # Save full AnnData with embedding in obsm
    h5ad_out = output_dir / f"{dataset_name}_scgpt_embeddings.h5ad"
    embedded.write(h5ad_out, compression="gzip")
    print(f"  - h5ad saved -> {h5ad_out}")

    # Save embedding matrix as CSV with metadata columns
    emb_matrix = embedded.obsm["X_scGPT"]
    df = pd.DataFrame(emb_matrix, columns=[f"dim_{i}" for i in range(emb_matrix.shape[1])])

    for col_idx, col in enumerate(config["meta_cols"]):
        if col in embedded.obs.columns:
            df.insert(col_idx, col, embedded.obs[col].values)

    csv_out = output_dir / f"{dataset_name}_scgpt_embeddings.csv"
    df.to_csv(csv_out, index=False)
    print(f"  - CSV saved  -> {csv_out}")


def main():
    OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

    for config in DATASETS_CONFIG:
        extract_embeddings(config)


if __name__ == "__main__":
    main()
