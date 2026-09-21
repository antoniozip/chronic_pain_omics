import os
import random

import numpy as np


def set_seed(seed: int = 42) -> None:
    """Set random seeds for reproducibility across Python, NumPy."""
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
