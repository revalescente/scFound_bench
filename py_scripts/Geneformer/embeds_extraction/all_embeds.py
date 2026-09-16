#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from datasets import load_from_disk
from transformers import BertForMaskedLM


def extract_all_hidden_states(config, model, device, output_base_dir, batch_size=64):
    dataset_name = config["name"]
    dataset_path = config["dataset_path"]
    meta_cols = config["meta_cols"]

    print(f"\n==========================================")
    print(f" Extracting ALL Layers for Dataset: {dataset_name.upper()}")
    print(f"==========================================")

    output_dir = Path(output_base_dir) / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset = load_from_disk(dataset_path)
    num_layers = model.config.num_hidden_layers  # 12 per Geneformer-V2-104M

    layer_embeddings = {l: [] for l in range(num_layers + 1)}
    metadata_collected = {col: [] for col in meta_cols}

    print(f"Estrazione hidden states per {len(dataset)} cellule ({num_layers + 1} layer totali)...")

    with torch.no_grad():
        for i in range(0, len(dataset), batch_size):
            batch = dataset[i : i + batch_size]

            # Padding dinamico per gestire lunghezze variabili nel batch
            batch_tensors = [torch.tensor(x, dtype=torch.long) for x in batch["input_ids"]]
            input_ids = torch.nn.utils.rnn.pad_sequence(
                batch_tensors, batch_first=True, padding_value=0
            ).to(device)

            attention_mask = (input_ids != 0).long()

            # Forward pass con estrazione di tutte le rappresentazioni interne
            outputs = model(input_ids=input_ids, attention_mask=attention_mask)
            hidden_states = outputs.hidden_states

            mask_expanded = (
                attention_mask.unsqueeze(-1).expand_as(hidden_states[0]).float()
            )

            # Mean pooling escludendo i token di padding
            for l, state in enumerate(hidden_states):
                sum_embeddings = torch.sum(state * mask_expanded, dim=1)
                sum_mask = mask_expanded.sum(dim=1).clamp(min=1e-9)
                mean_pooled = (sum_embeddings / sum_mask).cpu().numpy()
                layer_embeddings[l].append(mean_pooled)

            # Estrazione dei metadati specifici
            for col in meta_cols:
                if col in batch:
                    metadata_collected[col].extend(batch[col])

    # Salvataggio di un file CSV distinto per ciascun layer
    for l in range(num_layers + 1):
        layer_mat = np.vstack(layer_embeddings[l])
        df = pd.DataFrame(layer_mat)

        # Inserimento colonne metadati all'inizio
        for col_idx, col in enumerate(meta_cols):
            if metadata_collected[col]:
                df.insert(col_idx, col, metadata_collected[col])

        out_file = output_dir / f"{dataset_name}_geneformer_layer_{l}.csv"
        df.to_csv(out_file, index=False)
        print(f"  - Layer {l} salvato -> {out_file}")


def main():
    model_dir = "/software/apps/huggingface/hub/models--ctheodoris--Geneformer/snapshots/04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5/Geneformer-V2-104M/"
    output_base_dir = "/projects/shared/intronic_bam/datasets/embeddings/geneformer_layers"

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Caricamento modello su {device}...")
    model = BertForMaskedLM.from_pretrained(
        model_dir, output_hidden_states=True
    ).to(device)
    model.eval()

    # Configurazione di tutti e 3 i dataset
    datasets_config = [
        {
            "name": "be1",
            "dataset_path": "/projects/shared/intronic_bam/datasets/geneformer/be1.dataset",
            "meta_cols": ["Barcode", "Sample"],
        },
        {
            "name": "cb",
            "dataset_path": "/projects/shared/intronic_bam/datasets/geneformer/cb.dataset",
            "meta_cols": ["Barcode", "celltype"],
        },
        {
            "name": "mix",
            "dataset_path": "/projects/shared/intronic_bam/datasets/geneformer/mix.dataset",
            "meta_cols": ["Barcode", "cell_line"],
        },
    ]

    # --- MODALITÀ SINGOLO DATASET (Esegue solo 'be1') ---
    extract_all_hidden_states(datasets_config[0], model, device, output_base_dir)

    # --- MODALITÀ MULTI-DATASET (Scommenta qui sotto per eseguire su TUTTI i dataset) ---
    # for config in datasets_config:
    #     extract_all_hidden_states(config, model, device, output_base_dir)


if __name__ == "__main__":
    main()
