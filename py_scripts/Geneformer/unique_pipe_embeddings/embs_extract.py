#!/usr/bin/env python3
from pathlib import Path
from geneformer import EmbExtractor


def extract_dataset_embeddings(config, model_dir, output_dir):
    dataset_name = config["name"]
    output_prefix = f"{dataset_name}_geneformer"
    input_data_dir = (
        f"/projects/shared/intronic_bam/datasets/geneformer/{dataset_name}.dataset"
    )

    print(f"\n==========================================")
    print(f" Extracting Embeddings: {dataset_name.upper()}")
    print(f"==========================================")

    embex = EmbExtractor(
        model_type="Pretrained",
        num_classes=config.get("num_classes", 0),
        emb_mode="cell",
        cell_emb_style="mean_pool",
        filter_data=None,
        max_ncells=None,
        emb_layer=-1,
        emb_label=config["emb_label"],
        labels_to_plot=config["labels_to_plot"],
        forward_batch_size=128,
        nproc=1,
        summary_stat=None,
    )

    embs_df = embex.extract_embs(
        model_directory=model_dir,
        input_data_file=input_data_dir,
        output_directory=output_dir,
        output_prefix=output_prefix,
        output_torch_embs=False,
    )

    print(f"--- Finished {dataset_name} ---")
    print(f"Embeddings saved to: {Path(output_dir) / f'{output_prefix}.csv'}")
    print(f"DataFrame shape: {embs_df.shape}")


def main():
    model_dir = "/software/apps/huggingface/hub/models--ctheodoris--Geneformer/snapshots/04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5/Geneformer-V2-104M/"
    output_dir = "/projects/shared/intronic_bam/datasets/embeddings/geneformer"

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Configurazione specifica per dataset
    datasets_config = [
        {
            "name": "be1",
            "num_classes": 0,
            "emb_label": ["Barcode", "Sample"],
            "labels_to_plot": ["Sample"],
        },
        {
            "name": "cb",
            "num_classes": 0,
            "emb_label": ["Barcode", "celltype"],
            "labels_to_plot": ["celltype"],
        },
        {
            "name": "mix",
            "num_classes": 0,
            "emb_label": ["Barcode", "cell_line"],
            "labels_to_plot": ["cell_line"],
        },
    ]

    for config in datasets_config:
        extract_dataset_embeddings(config, model_dir, output_dir)


if __name__ == "__main__":
    main()
