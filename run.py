import random
import numpy as np
import torch

from pipeline import run

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)
torch.cuda.manual_seed_all(42)

torch.backends.cudnn.benchmark = False
torch.use_deterministic_algorithms(True)

run("geneformer", ["be1", "sc_mixology", "cord_blood"], batch_size=32)
run("scgpt", ["be1", "sc_mixology", "cord_blood"], batch_size=32)