"""
analysis.types
==============
Shared type aliases used across the analysis package.
"""

from numpy.typing import NDArray
import numpy as np
from dataclasses import dataclass,field

from typing import Dict, Optional, Tuple

# Canonical position dict type: droplet_id -> (N, 2) float array
PositionsDict = Dict[int, NDArray[np.floating]]

@dataclass(slots=True)
class PairDistanceData:
    distances: NDArray[np.float64]
    pair_map: list[tuple[int, int]]


    _mean: Optional[NDArray[np.float64]] = field(default=None, init=False)
    _std: Optional[NDArray[np.float64]] = field(default=None, init=False)

    @property
    def npairs(self) -> int:
        return len(self.pair_map)

    @property
    def nframes(self) -> int:
        return self.distances.shape[0]
    
    def pair(self, k: int) -> tuple[int, int]:
        return self.pair_map[k]

    @property
    def mean(self):
        if self._mean is None:
            self._mean = self.distances.mean(axis=0)
        return self._mean


    @property
    def std(self)->NDArray[np.float64]:
        if self._std is None:
            self._std = self.distances.std(axis=0)
            
            # self._std[self._std == 0] = 1.0
            eps = np.finfo(np.float64).eps
            self._std = np.maximum(self._std, eps)
            
        return self._std

    def standardized(self) -> NDArray[np.float64]:
        return np.asarray(
            (self.distances - self.mean) / self.std,
            dtype=np.float64
        )
    
    def destandardize(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Convert standardized pair distances back to physical units.
        """
        return x * self.std + self.mean

    def distance_change(self, dz: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Convert standardized distance changes into physical distance changes.
        """
        return dz * self.std