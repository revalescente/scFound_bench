import anndata as ad
import numpy as np
from scipy import sparse
from torch.utils.data import Dataset

from config import DATASETS


class SingleCellDataset(Dataset):
    """Read QC-filtered raw counts from an AnnData file, one cell at a time."""

    def __init__(self, name, split="all"):
        cfg = DATASETS[name]
        self.adata = ad.read_h5ad(cfg["path"], backed="r")
        obs = self.adata.obs
        mask = ~obs["discard"].astype(bool).to_numpy()
        if split != "all":
            train = obs["is_train"].astype(str).str.upper().to_numpy() == "TRUE"
            mask &= train if split == "train" else ~train

        self.indices = np.flatnonzero(mask)
        self.labels = obs[cfg["label"]].iloc[self.indices].tolist()
        self.cell_ids = obs.index[self.indices].tolist()
        var = self.adata.var
        self.gene_names = var[cfg["genes"]].astype(str).tolist() if cfg["genes"] else var.index.tolist()
        self.ensembl_ids = var[cfg["ensembl"]].astype(str).tolist() if cfg["ensembl"] else None

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        row = self.adata.layers["counts"][int(self.indices[idx]), :]
        return row.toarray().ravel() if sparse.issparse(row) else np.asarray(row).ravel()

    def close(self):
        self.adata.file.close()
