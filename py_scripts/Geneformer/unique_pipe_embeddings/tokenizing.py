#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import scanpy as sc
from geneformer import TranscriptomeTokenizer
import mygene


def fetch_ensembl_ids(gene_symbols):
    """Mappa una lista di Simboli genici in Ensembl ID usando MyGeneInfo."""
    mg = mygene.MyGeneInfo()
    print("  - Scaricamento degli Ensembl ID via MyGene...")
    results = mg.querymany(
        gene_symbols, scopes="symbol", fields="ensembl.gene", species="human"
    )

    ensembl_map = {}
    for item in results:
        query_sym = item["query"]
        if query_sym in ensembl_map:
            continue

        if "ensembl" in item:
            if isinstance(item["ensembl"], list):
                ensembl_map[query_sym] = item["ensembl"][0]["gene"]
            else:
                ensembl_map[query_sym] = item["ensembl"]["gene"]
        else:
            ensembl_map[query_sym] = None

    return ensembl_map


def process_and_tokenize(config, base_output_dir):
    prefix = config["prefix"]
    input_h5ad = config["input_h5ad"]
    prep_dir = Path(config["prep_dir"])

    prep_dir.mkdir(parents=True, exist_ok=True)
    prepped_h5ad = prep_dir / f"{prefix}_prepped.h5ad"

    print(f"\n==========================================")
    print(f" Processing Dataset: {prefix}")
    print(f"==========================================")
    print(f"Loading dataset: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)

    print("Formatting metadata for Geneformer...")

    # 1. Assegnazione/Mappatura Ensembl ID
    if config["ensembl_source"] == "var_col":
        adata.var["ensembl_id"] = adata.var[config["var_col_name"]]
    elif config["ensembl_source"] == "mygene":
        ensembl_map = fetch_ensembl_ids(adata.var.index.tolist())
        adata.var["ensembl_id"] = adata.var.index.map(ensembl_map)

    # 2. Assegnazione n_counts
    counts_src = config["counts_source"]
    if counts_src == "calc_sum":
        counts = adata.X.sum(axis=1)
        adata.obs["n_counts"] = np.asarray(counts).flatten()
    else:
        adata.obs["n_counts"] = np.asarray(adata.obs[counts_src]).flatten()

    # 3. Gestione colonna Barcode
    if "Barcode" not in adata.obs:
        adata.obs["Barcode"] = adata.obs.index

    print(f"Saving prepped dataset to: {prepped_h5ad}")
    adata.write_h5ad(prepped_h5ad)

    print("Tokenizing data...")
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict=config["custom_attr_dict"],
        nproc=10,
    )

    tk.tokenize_data(
        data_directory=str(prep_dir),
        output_directory=str(base_output_dir),
        output_prefix=prefix,
        file_format="h5ad",
        use_generator=False,
    )

    tokenized_path = base_output_dir / f"{prefix}.dataset"
    print(f"Dataset {prefix} tokenized successfully at: {tokenized_path}")


def main():
    output_dir = Path("/projects/shared/intronic_bam/datasets/geneformer")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Configurazione per ciascun dataset
    datasets_config = [
        {
            "prefix": "be1",
            "input_h5ad": "/projects/shared/intronic_bam/datasets/anndata/be1_scGPT_embeddings.h5ad",
            "prep_dir": "/projects/shared/intronic_bam/datasets/geneformer/prep_be1",
            "ensembl_source": "var_col",
            "var_col_name": "ID",
            "counts_source": "calc_sum",
            "custom_attr_dict": {"Sample": "Sample", "Barcode": "Barcode"},
        },
        {
            "prefix": "cb",
            "input_h5ad": "/projects/shared/intronic_bam/datasets/anndata/cb_scGPT_embeddings.h5ad",
            "prep_dir": "/projects/shared/intronic_bam/datasets/geneformer/prep_cb",
            "ensembl_source": "mygene",
            "counts_source": "total",
            "custom_attr_dict": {"celltype": "celltype", "Barcode": "Barcode"},
        },
        {
            "prefix": "mix",
            "input_h5ad": "/projects/shared/intronic_bam/datasets/anndata/mix_scGPT_embeddings.h5ad",
            "prep_dir": "/projects/shared/intronic_bam/datasets/geneformer/prep_mix",
            "ensembl_source": "mygene",
            "counts_source": "total_counts",
            "custom_attr_dict": {"cell_line": "cell_line", "Barcode": "Barcode"},
        },
    ]

    for config in datasets_config:
        process_and_tokenize(config, output_dir)


if __name__ == "__main__":
    main()
