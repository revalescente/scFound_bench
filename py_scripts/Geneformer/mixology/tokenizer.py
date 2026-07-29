#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import scanpy as sc
from geneformer import TranscriptomeTokenizer
import mygene


def main():
    input_h5ad = "/projects/shared/intronic_bam/datasets/anndata/mix_scGPT_embeddings.h5ad"
    prep_dir = "/projects/shared/intronic_bam/datasets/geneformer/prep_mix"
    output_dir = "/projects/shared/intronic_bam/datasets/geneformer"
    output_prefix = "mix"

    # Assicurati che le directory esistano
    Path(prep_dir).mkdir(parents=True, exist_ok=True)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    prepped_h5ad = Path(prep_dir) / "mix_prepped.h5ad"

    print(f"Loading original dataset: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)

    print("Formatting metadata for Geneformer...")
    # 1. Procuro gli ensembl_id dal client
    # 1a. Inizializza il client MyGene
    mg = mygene.MyGeneInfo()

    # 1b. Estrai la lista di tutti i simboli dei geni dall'indice di adata.var
    gene_list = adata.var.index.tolist()

    # 1c. Esegui la query in blocco per tutti i geni
    print("Scaricamento degli Ensembl ID in corso...")
    results = mg.querymany(gene_list, scopes='symbol', fields='ensembl.gene', species='human')

    # 1d. Crea un dizionario per mappare il Simbolo all'Ensembl ID in modo sicuro
    ensembl_map = {}
    for item in results:
        query_sym = item['query']

        # Se abbiamo già mappato questo gene (in caso di risultati multipli), teniamo il primo
        if query_sym in ensembl_map:
            continue

        # Se la API ha trovato una corrispondenza con Ensembl
        if 'ensembl' in item:
            # A volte il risultato è una lista (se un simbolo mappa su più Ensembl ID)
            if isinstance(item['ensembl'], list):
                ensembl_map[query_sym] = item['ensembl'][0]['gene']
            else:
                ensembl_map[query_sym] = item['ensembl']['gene']
        else:
            # Se non esiste l'ID per quel gene, mettiamo None (diventerà NaN in Pandas)
            ensembl_map[query_sym] = None

    # 1e. Mappa il dizionario all'indice del tuo AnnData per creare la nuova colonna
    adata.var['ensembl_id'] = adata.var.index.map(ensembl_map)


    # 2. Calcola i conteggi totali per cellula (conversione sicura in array 1D)
    adata.obs["n_counts"] = adata.obs["total_counts"]

    #3 Inserimento colonna Barcode
    adata.obs["Barcode"] = adata.obs.index

    print(f"Saving prepped dataset to {prepped_h5ad}")
    adata.write_h5ad(prepped_h5ad)

    print("Initializing TranscriptomeTokenizer...")
    # Preserviamo sia Sample che Barcode nei metadati del dataset tokenizzato
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"cell_line": "cell_line", "Barcode": "Barcode"},
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
