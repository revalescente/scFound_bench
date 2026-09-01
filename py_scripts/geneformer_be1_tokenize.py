import logging
import scanpy as sc
import anndata as ad
from pathlib import Path
from geneformer import TranscriptomeTokenizer

logger = logging.getLogger(__name__)

def main():
    # Define paths
    input_h5ad = "/projects/shared/intronic_bam/datasets/anndata/be1.h5ad"
    # Directory to store the intermediate correctly formatted h5ad file
    prep_dir = "/projects/shared/intronic_bam/datasets/geneformer/prep_be1"
    output_dir = "/projects/shared/intronic_bam/datasets/geneformer"
    output_prefix = "be1"
    
    Path(prep_dir).mkdir(parents=True, exist_ok=True)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    prepped_h5ad = Path(prep_dir) / "be1_prepped.h5ad"
    
    print(f"Loading original dataset: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    
    print("Formatting metadata for Geneformer...")
    # Geneformer TranscriptomeTokenizer hardcodes the requirement for:
    # 1. 'ensembl_id' in adata.var
    # 2. 'n_counts' in adata.obs
    
    # Assuming 'ID' contains the Ensembl IDs. If it's 'Symbol' or another column, adjust here.
    if 'ID' in adata.var.columns:
        adata.var['ensembl_id'] = adata.var['ID']
    else:
        # Fallback if 'ID' isn't available, maybe gene names/index are Ensembl IDs
        adata.var['ensembl_id'] = adata.var.index

    # Assuming 'sum' or 'total' contains the total read counts per cell
    if 'sum' in adata.obs.columns:
        adata.obs['n_counts'] = adata.obs['sum']
    elif 'total' in adata.obs.columns:
        adata.obs['n_counts'] = adata.obs['total']
    else:
        # If neither exists, calculate it from the raw counts matrix
        adata.obs['n_counts'] = adata.X.sum(axis=1)
        
    print(f"Saving prepped dataset to {prepped_h5ad}")
    # Save the prepared anndata
    adata.write_h5ad(prepped_h5ad)
    
    # Now we can tokenize
    print("Initializing TranscriptomeTokenizer...")
    # Add custom attributes we want to retain in the tokenized dataset.
    # We map output_name: input_name
    # From earlier notebooks, 'Sample' is the cell type key for BE1.
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"Sample": "Sample"}, 
        nproc=16
    )
    
    print("Tokenizing data...")
    tk.tokenize_data(
        data_directory=prep_dir,
        output_directory=output_dir,
        output_prefix=output_prefix,
        file_format="h5ad",
        use_generator=False
    )
    
    print(f"Tokenization complete! Tokenized dataset saved at: {Path(output_dir) / (output_prefix + '.dataset')}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
