"""
analysis.types
==============
Shared type aliases used across the analysis package.
"""

from numpy.typing import NDArray
import numpy as np
from typing import Dict, Optional, Tuple

# Canonical position dict type: droplet_id -> (N, 2) float array
PositionsDict = Dict[int, NDArray[np.floating]]
