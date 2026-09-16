import argparse
import glob
import os
from concurrent.futures import ProcessPoolExecutor

import anndata as ad


def parse_args():
    parser = argparse.ArgumentParser(description="Splitta le partizioni CellxGene per dataset_id")
    parser.add_argument("--tissue-dir", type=str, required=True, help="Directory del tessuto (es. /projects/shared/cellxgene/blood)")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory di output finale per i dataset")
    parser.add_argument("--n-workers", type=int, default=8, help="Numero di processori paralleli")
    return parser.parse_args()

def process_single_dataset(dataset_id, file_slices, output_dir):
    out_file = os.path.join(output_dir, f"{dataset_id}.h5ad")
    tmp_file = os.path.join(output_dir, f"{dataset_id}.h5ad.tmp")

    # 1. Calcola il numero totale di cellule attese per questo dataset_id
    expected_n_obs = sum(len(indices) for _, indices in file_slices)

    # 2. CHECK INTEGRITÀ: Se il file esiste già, verifichiamo se è completo
    if os.path.exists(out_file):
        try:
            existing_adata = ad.read_h5ad(out_file, backed="r")
            existing_n_obs = existing_adata.n_obs
            existing_adata.file.close()

            if existing_n_obs == expected_n_obs:
                print(f"[SKIP] {dataset_id} già completo ({expected_n_obs} cellule).")
                return
            else:
                print(f"[RETRY] {dataset_id} incompleto! Attese {expected_n_obs}, trovate {existing_n_obs}. Rigenero...")
        except Exception as e:
            print(f"[RETRY] {dataset_id} corrotto ({e}). Rigenero...")

    # 3. Estrazione dati dalle partizioni caricando la slice in RAM
    adatas = []
    for fpath, cell_indices in file_slices:
        adata_backed = ad.read_h5ad(fpath, backed="r")
        sub_adata = adata_backed[cell_indices].to_memory()
        adata_backed.file.close()
        adatas.append(sub_adata)

    # Concatena le frazioni
    if len(adatas) == 1:
        merged_adata = adatas[0]
    else:
        merged_adata = ad.concat(adatas, merge="same")

    # Se ci sono indici duplicati nelle cellule (.obs), li rende unici
    if not merged_adata.obs.index.is_unique:
        merged_adata.obs_names_make_unique()
    merged_adata.obs.index.name = None  # Previeni conflitti con colonne esistenti

    # Se ci sono indici duplicati nei geni (.var), li rende unici
    if not merged_adata.var.index.is_unique:
        merged_adata.var_names_make_unique()
    merged_adata.var.index.name = None

    # 5. SCRITTURA ATOMICA (.tmp -> .h5ad)
    merged_adata.write_h5ad(tmp_file, compression="gzip")
    os.replace(tmp_file, out_file)
    print(f"[DONE] Salvato {dataset_id} ({merged_adata.n_obs} cellule) -> {out_file}")

def main():
    args = parse_args()
    tissue_name = os.path.basename(os.path.normpath(args.tissue_dir))
    out_tissue_dir = os.path.join(args.output_dir, tissue_name)
    os.makedirs(out_tissue_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.tissue_dir, "*.h5ad")))
    print(f"=== Elaborazione Tessuto: {tissue_name} ({len(files)} partizioni) ===")

    # FASE 1: Scan rapido dei metadati (.obs)
    print("Fase 1: Mappatura dei dataset_id nelle partizioni...")
    dataset_map = {} # { dataset_id: [(file_path, [cell_indices]), ... ] }

    for fpath in files:
        # backed='r' carica SOLTANTO i metadati in frazioni di secondo
        adata_backed = ad.read_h5ad(fpath, backed="r")
        dataset_series = adata_backed.obs["dataset_id"]

        # Raggruppa gli indici per ciascun dataset_id presente in questa partizione
        grouped = dataset_series.groupby(dataset_series, observed=True).groups
        for ds_id, idx_label in grouped.items():
            # Converte le etichette di indice in posizioni numeriche di riga
            pos_indices = adata_backed.obs.index.get_indexer(idx_label)
            if ds_id not in dataset_map:
                dataset_map[ds_id] = []
            dataset_map[ds_id].append((fpath, pos_indices))

        adata_backed.file.close()

    print(f"Trovati {len(dataset_map)} dataset_id unici in {tissue_name}.")

    # FASE 2: Estrazione e Scrittura Parallela su CPU
    print(f"Fase 2: Estrazione e salvataggio in parallelo con {args.n_workers} worker...")

    with ProcessPoolExecutor(max_workers=args.n_workers) as executor:
        futures = [
            executor.submit(process_single_dataset, ds_id, file_slices, out_tissue_dir)
            for ds_id, file_slices in dataset_map.items()
        ]
        for future in futures:
            future.result()

if __name__ == "__main__":
    main()
