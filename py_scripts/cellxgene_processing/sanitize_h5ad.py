import os
import glob
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import anndata as ad

def parse_args():
    parser = argparse.ArgumentParser(description="Bonifica completa dell'intero database .h5ad per compatibilità con R")
    parser.add_argument(
        "--base-dir",
        type=str,
        default="/projects/shared/cellxgene_split_by_dataset",
        help="Directory radice da scansionare"
    )
    parser.add_argument("--n-workers", type=int, default=16, help="Numero di worker paralleli")
    return parser.parse_args()

def check_and_fix_file(fpath):
    """
    Verifica se un file .h5ad ha indici duplicati in .obs o .var.
    Se li trova, rende unici gli indici senza creare conflitti con i nomi delle colonne.
    """
    rel_path = os.path.relpath(fpath)

    try:
        # FASE 1: Check rapido senza caricare la matrice in RAM
        adata_check = ad.read_h5ad(fpath, backed="r")
        obs_unique = adata_check.obs.index.is_unique
        var_unique = adata_check.var.index.is_unique
        adata_check.file.close()

        # Se sono già unici, saltiamo la modifica
        if obs_unique and var_unique:
            return ("OK", rel_path)

        # FASE 2: Correzione del file con problemi
        print(f"[FIX] Trovati indici duplicati in {rel_path}. Correzione in corso...")

        adata = ad.read_h5ad(fpath) # Caricamento in RAM

        # 1. Correzione indici cellule (.obs)
        if not obs_unique:
            adata.obs_names_make_unique()
            adata.obs.index.name = None  # Evita conflitti tra nome indice e colonne

        # 2. Correzione indici geni (.var)
        if not var_unique:
            adata.var_names_make_unique()
            adata.var.index.name = None

        # Scrittura atomica via file temporaneo
        tmp_file = f"{fpath}.sanitize.tmp"
        adata.write_h5ad(tmp_file, compression="gzip")
        os.replace(tmp_file, fpath)

        return ("FIXED", rel_path)

    except Exception as e:
        return ("ERROR", f"{rel_path}: {str(e)}")

def main():
    args = parse_args()

    print(f"=== AVVIO BONIFICA DATABASE: {args.base_dir} ===")

    search_pattern = os.path.join(args.base_dir, "**", "*.h5ad")
    all_files = sorted(glob.glob(search_pattern, recursive=True))

    if not all_files:
        print(f"Nessun file .h5ad trovato in {args.base_dir}")
        return

    print(f"Trovati {len(all_files)} file .h5ad da verificare.")
    print(f"Esecuzione in parallelo con {args.n_workers} worker...\n")

    stats = {"OK": 0, "FIXED": 0, "ERROR": 0}

    with ProcessPoolExecutor(max_workers=args.n_workers) as executor:
        futures = {executor.submit(check_and_fix_file, fpath): fpath for fpath in all_files}

        for future in as_completed(futures):
            status, msg = future.result()
            stats[status] += 1

            if status == "FIXED":
                print(f"  ✅ CORRETTO: {msg}")
            elif status == "ERROR":
                print(f"  ❌ ERRORE: {msg}")

    print("\n" + "=" * 50)
    print("=== RIEPILOGO BONIFICA ===")
    print(f"File già corretti (OK) : {stats['OK']}")
    print(f"File modificati (FIXED): {stats['FIXED']}")
    print(f"File con errore (ERROR): {stats['ERROR']}")
    print("=" * 50)

if __name__ == "__main__":
    main()
