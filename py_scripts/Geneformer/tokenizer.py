#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import scanpy as sc
from geneformer import TranscriptomeTokenizer


def main():
    input_h5ad = "/projects/shared/intronic_bam/datasets/anndata/be1_scGPT_embeddings.h5ad"
    prep_dir = "/projects/shared/intronic_bam/datasets/geneformer/prep_be1"
    output_dir = "/projects/shared/intronic_bam/datasets/geneformer"
    output_prefix = "be1"

    # Assicurati che le directory esistano
    Path(prep_dir).mkdir(parents=True, exist_ok=True)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    prepped_h5ad = Path(prep_dir) / "be1_prepped.h5ad"

    print(f"Loading original dataset: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)

    print("Formatting metadata for Geneformer...")
    # 1. Assegna gli Ensembl IDs alla colonna richiesta da Geneformer
    adata.var["ensembl_id"] = adata.var["ID"]

    # 2. Calcola i conteggi totali per cellula (conversione sicura in array 1D)
    counts = adata.X.sum(axis=1)
    adata.obs["n_counts"] = np.asarray(counts).flatten()

    print(f"Saving prepped dataset to {prepped_h5ad}")
    adata.write_h5ad(prepped_h5ad)

    print("Initializing TranscriptomeTokenizer...")
    # Preserviamo sia Sample che Barcode nei metadati del dataset tokenizzato
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"Sample": "Sample", "Barcode": "Barcode"},
        nproc=10,
    )

    print("Tokenizing data...")
    tk.tokenize_data(
        data_directory=prep_dir,
        output_directory=output_dir,
        output_prefix=output_prefix,
        file_format="h5ad",
        use_generator=False,
    )

    tokenized_path = Path(output_dir) / f"{output_prefix}.dataset"
    print(f"Tokenization complete! Saved at: {tokenized_path}")


if __name__ == "__main__":
    main()
