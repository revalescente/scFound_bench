import glob
import os
import anndata as ad

# Percorso radice contenente tutte le cartelle dei dataset
base_dir = "/projects/shared/cellxgene"

# Trova tutte le sottocartelle escludendo 'faiss_index'
subdirs = sorted([
    d for d in os.listdir(base_dir)
    if os.path.isdir(os.path.join(base_dir, d)) and d != "faiss_index"
])

print(f"Trovati {len(subdirs)} dataset da analizzare: {subdirs}\n")

for subdir in subdirs:
    data_dir = os.path.join(base_dir, subdir)
    files = sorted(glob.glob(os.path.join(data_dir, "*.h5ad")))

    print("#" * 80)
    print(f" DATASET: {subdir.upper()} ({len(files)} file .h5ad trovati)")
    print("#" * 80 + "\n")

    if not files:
        print(f"Nessun file .h5ad trovato nella cartella '{data_dir}'\n")
        continue

    for fpath in files:
        fname = os.path.basename(fpath)
        print("=" * 70)
        print(f"FILE: {subdir}/{fname}")
        print("=" * 70)

        try:
            # Legge in modalità backed ('r') per caricare solo i metadati
            adata = ad.read_h5ad(fpath, backed="r")

            # 1. Dimensioni
            n_obs, n_vars = adata.shape
            print(f"Dimensioni: {n_obs} cellule (obs) x {n_vars} geni (var)")

            # 2. Nomi variabili in .obs e .var
            print("\n--- .obs (variabili cellule) ---")
            print(list(adata.obs.columns))

            print("\n--- .var (variabili geni) ---")
            print(list(adata.var.columns))

            # 3. Altri oggetti presenti
            print("\n--- Altri oggetti ---")
            slots = {
                "obsm": list(adata.obsm.keys()),
                "varm": list(adata.varm.keys()),
                "uns": list(adata.uns.keys()),
                "obsp": list(adata.obsp.keys()),
                "varp": list(adata.varp.keys()),
                "layers": list(adata.layers.keys()),
            }

            for slot_name, keys in slots.items():
                if keys:
                    print(f"  - {slot_name}: {keys}")
                else:
                    print(f"  - {slot_name}: [Vuoto]")

            # Chiude il file aperto in modalità backed
            if adata.isbacked:
                adata.file.close()

        except Exception as e:
            print(f"Errore nella lettura del file {fname}: {e}")

        print("\n")
