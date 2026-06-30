"""
analysis.types
==============
Shared type aliases used across the analysis package.
"""

from numpy.typing import NDArray
import numpy as np
from dataclasses import dataclass

from typing import Dict, Optional, Tuple

# Canonical position dict type: droplet_id -> (N, 2) float array
PositionsDict = Dict[int, NDArray[np.floating]]

@dataclass(slots=True)
class PairDistanceData:
    distances: NDArray[np.float64]
    pair_map: list[tuple[int, int]]

    @property
    def npairs(self) -> int:
        return len(self.pair_map)

    @property
    def nframes(self) -> int:
        return self.distances.shape[0]
    
    def pair(self, k: int) -> tuple[int, int]:
        return self.pair_map[k]
    
    def standardized(self) -> NDArray[np.float64]:
        mean = self.distances.mean(axis=0)
        std = self.distances.std(axis=0)

        return (self.distances - mean) / std