from pathlib import Path

DATASETS = {
    "be1": {
        "path": "/projects/shared/intronic_bam/datasets/anndata/be1.h5ad",
        "label": "Sample", "genes": "Symbol", "ensembl": "ID",
    },
    "sc_mixology": {
        "path": "/projects/shared/intronic_bam/datasets/anndata/sc_mix.h5ad",
        "label": "cell_line", "genes": None, "ensembl": None,
    },
    "cord_blood": {
        "path": "/projects/shared/intronic_bam/datasets/anndata/cord_blood.h5ad",
        "label": "celltype", "genes": "gene_names", "ensembl": None,
    },
}

SCGPT_DIR = Path(
    "/projects/shared/intronic_bam/scGPT_supplements/save/scGPT_human"
)
GENEFORMER_DIR = Path(
    "/software/apps/huggingface/hub/models--ctheodoris--Geneformer/"
    "snapshots/04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5/Geneformer-V2-104M"
)
GENEFORMER_ASSETS = Path(
    "/software/apps/huggingface/hub/models--ctheodoris--Geneformer/"
    "snapshots/1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5/geneformer"
)

OUTPUT_DIR = Path("outputs3")
