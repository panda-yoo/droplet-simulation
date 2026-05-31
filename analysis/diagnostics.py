"""
analysis.diagnostics
====================
Ensemble bias diagnostics.

These functions quantify the effect of the track-length filter that retains
only full-length trajectories for ensemble averaging. They do NOT alter the
analysis pipeline; they are purely informational.

Contents
--------
- retained_vs_discarded_lengths    : track-length distributions
- retained_vs_discarded_speed_stats: speed statistics for retained vs discarded
- pair_weighted_msd                : alias → analysis.msd.compute_pair_weighted_msd
- pair_weighted_vacf               : pair-weighted VACF across all tracks
"""

import numpy as np
from numpy.typing import NDArray
from typing import Dict, List, Tuple

from analysis.types import PositionsDict
from analysis.msd import compute_pair_weighted_msd
from analysis.vacf import compute_velocity


# =========================================================
# Re-export pair_weighted_msd under the diagnostics namespace
# =========================================================

# Callers may import either from analysis.msd or analysis.diagnostics.
pair_weighted_msd = compute_pair_weighted_msd


# =========================================================
# Track length distributions
# =========================================================

def retained_vs_discarded_lengths(
    lengths: Dict[int, int],
) -> Tuple[List[int], List[int], int, int]:
    """
    Split track lengths into retained (full-length) and discarded (partial) groups.

    Parameters
    ----------
    lengths : dict
        Mapping droplet_id → track length (number of frames).
        As returned by ``load_positions_with_lengths``.

    Returns
    -------
    retained_lengths : list of int
        Lengths of full-length tracks (equal to the global maximum length).
    discarded_lengths : list of int
        Lengths of partial tracks (shorter than the global maximum).
    n_retained : int
        Count of retained tracks.
    n_discarded : int
        Count of discarded tracks.

    Notes
    -----
    A track is considered full-length if its length equals
    ``max(lengths.values())``.
    """
    if not lengths:
        return [], [], 0, 0

    max_len = max(lengths.values())
    retained = [l for l in lengths.values() if l == max_len]
    discarded = [l for l in lengths.values() if l != max_len]

    return retained, discarded, len(retained), len(discarded)


# =========================================================
# Speed statistics: retained vs discarded
# =========================================================

def retained_vs_discarded_speed_stats(
    positions_all: PositionsDict,
    lengths: Dict[int, int],
    dt: float = 1.0,
) -> Dict[str, float]:
    """
    Compute mean and standard deviation of time-averaged speed for retained
    and discarded tracks separately.

    Parameters
    ----------
    positions_all : PositionsDict
        All trajectories (including partial tracks), as returned by
        ``load_positions_with_lengths``.
    lengths : dict
        Track length per droplet ID.
    dt : float
        Time step between consecutive frames.

    Returns
    -------
    stats : dict with keys
        - ``retained_mean_speed``   : float
        - ``retained_std_speed``    : float
        - ``discarded_mean_speed``  : float
        - ``discarded_std_speed``   : float
        - ``n_retained``            : int
        - ``n_discarded``           : int

    Notes
    -----
    Speed is the mean over time of ||v(t)|| for each individual track,
    computed from ``compute_velocity`` (forward difference).
    Returns NaN for groups that contain no tracks.
    """
    if not lengths:
        return {
            "retained_mean_speed": np.nan,
            "retained_std_speed": np.nan,
            "discarded_mean_speed": np.nan,
            "discarded_std_speed": np.nan,
            "n_retained": 0,
            "n_discarded": 0,
        }

    max_len = max(lengths.values())
    retained_speeds: List[float] = []
    discarded_speeds: List[float] = []

    for did, r in positions_all.items():
        v = compute_velocity(r, dt)
        mean_speed = float(np.mean(np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2)))
        if lengths[did] == max_len:
            retained_speeds.append(mean_speed)
        else:
            discarded_speeds.append(mean_speed)

    def _safe_mean(lst: List[float]) -> float:
        return float(np.mean(lst)) if lst else np.nan

    def _safe_std(lst: List[float]) -> float:
        return float(np.std(lst, ddof=1)) if len(lst) > 1 else np.nan

    return {
        "retained_mean_speed": _safe_mean(retained_speeds),
        "retained_std_speed": _safe_std(retained_speeds),
        "discarded_mean_speed": _safe_mean(discarded_speeds),
        "discarded_std_speed": _safe_std(discarded_speeds),
        "n_retained": len(retained_speeds),
        "n_discarded": len(discarded_speeds),
    }


# =========================================================
# Pair-weighted VACF
# =========================================================

def pair_weighted_vacf(
    positions_dict: PositionsDict,
    dt: float = 1.0,
    lag_fraction: float = 0.3,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]:
    """
    Compute pair-weighted VACF across all trajectories (including partial tracks).

    Parameters
    ----------
    positions_dict : PositionsDict
        Mapping droplet_id → (N_i, 2) position array.
        Trajectories need not have equal length.
    dt : float
        Time step between consecutive frames.
    lag_fraction : float
        Fraction of the shortest trajectory length used as the maximum lag.

    Returns
    -------
    tau : ndarray, shape (L,)
        Lag times in frame units: 0, 1, …, L-1.
    vacf_weighted : ndarray, shape (L,)
        Pair-weighted VACF: total velocity dot-product / total pair count
        at each lag.

    Notes
    -----
    Unlike ``ensemble_vacf_analysis``, which averages per-trajectory VACFs,
    this estimator weights each velocity pair equally. It is intended for
    bias diagnosis only; use ``ensemble_vacf_analysis`` for all physics
    results.
    """
    if len(positions_dict) == 0:
        return np.array([], dtype=int), np.array([], dtype=float)

    lengths = [len(r) for r in positions_dict.values()]
    min_len = min(lengths)
    max_lag = max(1, int(lag_fraction * (min_len - 1)))

    corr_sum = np.zeros(max_lag + 1, dtype=float)
    pair_counts = np.zeros(max_lag + 1, dtype=float)

    for r in positions_dict.values():
        v = compute_velocity(r, dt)
        nv = len(v)
        for tau in range(min(max_lag + 1, nv)):
            dot = v[:nv - tau, 0] * v[tau:, 0] + v[:nv - tau, 1] * v[tau:, 1]
            corr_sum[tau] += np.sum(dot)
            pair_counts[tau] += len(dot)

    with np.errstate(invalid="ignore", divide="ignore"):
        vacf_weighted = corr_sum / pair_counts

    tau = np.arange(max_lag + 1)
    return tau, vacf_weighted
