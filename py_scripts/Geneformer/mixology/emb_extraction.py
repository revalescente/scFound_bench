#!/usr/bin/env python3
from geneformer import EmbExtractor
from pathlib import Path



def main():
    # Define paths for model, input dataset, and output directory
    model_dir = "/software/apps/huggingface/hub/models--ctheodoris--Geneformer/snapshots/04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5/Geneformer-V2-104M/"
    input_data_dir = (
        "/projects/shared/intronic_bam/datasets/geneformer/mix.dataset"
    )
    output_dir = "/projects/shared/intronic_bam/datasets/embeddings/geneformer"
    output_prefix = "mix_geneformer"

    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Initialize the EmbExtractor
    print("Initializing EmbExtractor...")
    embex = EmbExtractor(
        model_type="Pretrained",  # Use "CellClassifier" if fine-tuned
        num_classes=0,  # Number of classes (if classifier)
        emb_mode="cell",
        cell_emb_style="mean_pool",
        filter_data=None,  # Filtering dict if needed, e.g., {"Sample": ["A549"]}
        max_ncells=None,  # None for all cells, or integer to downsample
        emb_layer=-1,  # -1 extracts 2nd to last layer (recommended)
        emb_label=[
            "Barcode",
            "cell_line",
        ],  # Preserves Barcode and Sample in the output dataframe
        labels_to_plot=["cell_line"],  # Label used to color UMAP/Heatmap
        forward_batch_size=128,  # Adjust based on GPU memory
        nproc=1,  # Number of CPU processes
        summary_stat=None,  # None outputs full embeddings per cell
    )

    # Extract the embeddings
    print("Extracting Geneformer embeddings for BE1 dataset...")
    embs_df = embex.extract_embs(
        model_directory=model_dir,
        input_data_file=input_data_dir,
        output_directory=output_dir,
        output_prefix=output_prefix,
        output_torch_embs=False,
    )

    print("\n--- Extraction Finished ---")
    print(f"Embeddings saved to: {output_dir}/{output_prefix}.csv")
    print(f"Embeddings DataFrame shape: {embs_df.shape}")


if __name__ == "__main__":
    main()
