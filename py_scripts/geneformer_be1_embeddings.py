import logging
from pathlib import Path
from geneformer import EmbExtractor

logger = logging.getLogger(__name__)

def main():
    # Define paths for the model, tokenized input data, and output directory
    # Note: Geneformer expects the input data to be pre-tokenized into a .dataset folder.
    # Replace these paths with the actual paths to your Geneformer model and tokenized BE1 dataset.
    model_dir = "/path/to/Geneformer/model"
    input_data_dir = "/projects/shared/intronic_bam/datasets/geneformer/be1.dataset"
    output_dir = "/projects/shared/intronic_bam/datasets/embeddings/geneformer"
    output_prefix = "be1_geneformer"
    
    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Initialize the EmbExtractor
    # The BE1 dataset uses "Sample" for its cell type labels.
    embex = EmbExtractor(
        model_type="Pretrained",     # Use "CellClassifier" if you fine-tuned a model
        num_classes=0,               # Set to number of classes if using a classifier
        emb_mode="cell",
        cell_emb_style="mean_pool",
        filter_data=None,            # Add filtering dict if needed, e.g., {"Sample": ["some_sample"]}
        max_ncells=None,             # Set to a number (e.g., 1000) to downsample, or None for all cells
        emb_layer=-1,                # -1 extracts the 2nd to last layer (recommended for pretrained)
        emb_label=["Sample"],        # Label to append to output and use for plotting
        labels_to_plot=["Sample"],   # Label used to color UMAP and Heatmap
        forward_batch_size=64,       # Adjust based on your GPU memory
        nproc=16,                    # Number of CPU processes
        summary_stat=None            # Output full embeddings. Change to "mean"/"median" if memory constrained
    )

    # Extract the embeddings
    print("Extracting Geneformer embeddings for BE1 dataset...")
    embs_df = embex.extract_embs(
        model_directory=model_dir,
        input_data_file=input_data_dir,
        output_directory=output_dir,
        output_prefix=output_prefix,
        output_torch_embs=False
    )
    
    print(f"Embeddings extracted and saved to {output_dir}/{output_prefix}.csv")
    print(f"Embeddings DataFrame shape: {embs_df.shape}")

    # Plot UMAP and Heatmap of the embeddings colored by Sample
    print("Plotting embeddings...")
    embex.plot_embs(
        embs=embs_df, 
        plot_style="umap",
        output_directory=output_dir,
        output_prefix=output_prefix,
        max_ncells_to_plot=5000,     # Downsample for plotting if the dataset is large
        kwargs_dict={"palette": "Set1", "size": 50}
    )
    
    # To plot a heatmap, we can call plot_embs again with plot_style="heatmap"
    # embex.plot_embs(
    #    embs=embs_df, 
    #    plot_style="heatmap",
    #    output_directory=output_dir,
    #    output_prefix=output_prefix,
    # )

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
