# Single-cell embeddings: scGPT, Geneformer

## Run

```bash
python run.py
```

Or run only one model/dataset:

```python
from pipeline import run
run("geneformer", ["be1"], batch_size=8)
run("scgpt", ["be1", "sc_mixology"], batch_size=32)
```

Outputs are written to `outputs/<dataset>_<split>_<model>.h5ad`. Each output has `obs['label']` and one `obsm['X_<model>_layer_N']` per encoder layer. The order of cells matches the filtered AnnData dataset.

## Files

- `config.py`: dataset paths, checkpoint paths, output directory
- `datasets.py`: common AnnData reader (QC and train/test split)
- `models/scgpt.py`: scGPT-specific vocabulary, binning and CLS embeddings
- `models/geneformer.py`: Geneformer ranking, dynamic padding and layer-wise mean pooling
- `models/__init__.py`: model registry
- `pipeline.py`: shared inference + output writing
- `run.py`: experiment selection

## Adding another model

1. Add `models/new_model.py` with a class exposing `__init__(...)` and `embed(data, batch_size)`.
2. Register it in `models/__init__.py`.
3. Add its loading branch in `pipeline.load_model` (or refactor loading into a registry of small factory functions).

`embed()` returns a NumPy array of shape `(n_cells, n_layers, hidden_size)` so the shared pipeline needs no changes.
