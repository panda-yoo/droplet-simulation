"""
analysis.compute — COMPATIBILITY SHIM
======================================
All implementations have been moved to the canonical per-quantity modules.
This file re-exports every public name so that existing code and notebooks
that imported from ``analysis.compute`` continue to work unchanged.

Do NOT add new implementations here.
"""

# ── MSD ──────────────────────────────────────────────────
from analysis.msd import (
    compute_msd,
    fit_msd_power_law,
    ensemble_msd_analysis,
    compute_local_msd_exponent,
    compute_pair_weighted_msd,
)

# ── VACF ─────────────────────────────────────────────────
from analysis.vacf import (
    compute_velocity as cal_velocity,  # canonical name: compute_velocity
    compute_vacf,
    fit_vacf_exponential,
    ensemble_vacf_analysis,
)

# ── Orientation ───────────────────────────────────────────
from analysis.orientation import (
    compute_orientation_corr,
    fit_orientation_persistence,
    ensemble_orientation_corr_analysis,
)

# ── Speed / heading (kept for backward compat) ───────────
import numpy as np
from numpy.typing import NDArray
from typing import Dict, Tuple
from analysis.types import PositionsDict
from analysis.vacf import compute_velocity


def ensemble_speed_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray[np.int_], NDArray[np.floating], NDArray[np.floating]]:
    """Backward-compatible ensemble speed analysis (delegates to canonical velocity)."""
    ids = list(positions_dict.keys())
    min_len = min(len(compute_velocity(positions_dict[i], dt)) for i in ids)

    ensemble_speed = np.zeros(min_len)
    ensemble_v = np.zeros((min_len, 2))

    for did in ids:
        v = compute_velocity(positions_dict[did], dt)[:min_len]
        ensemble_speed += np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2)
        ensemble_v += v

    ensemble_speed /= len(ids)
    ensemble_v /= len(ids)
    t = np.arange(min_len)
    return t, ensemble_speed, ensemble_v


def orientation_analysis(
    ensemble_v: NDArray[np.floating],
) -> NDArray[np.floating]:
    """Backward-compatible heading angle from ensemble velocity vector."""
    theta = np.arctan2(ensemble_v[:, 1], ensemble_v[:, 0])
    return np.unwrap(theta)


# ── PSD ───────────────────────────────────────────────────
from analysis.psd import (
    compute_psd_components,
    fit_psd_powerlaw,
    ensemble_psd_analysis,
)