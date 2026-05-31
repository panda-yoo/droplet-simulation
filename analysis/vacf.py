"""
analysis.vacf
=============
Canonical velocity autocorrelation function (VACF) computations.

Velocity convention (used everywhere in the analysis package):
    v[i] = (r[i+1] - r[i]) / dt      (forward difference, shape (N-1, 2))

All time axes use frame-index integers internally.
Convert to seconds in callers with ``tau_seconds = tau * dt``.

Contents
--------
- compute_velocity      : forward-difference velocity from positions
- compute_vacf          : single-droplet VACF
- fit_vacf_exponential  : exponential fit → velocity persistence time τ_p
- ensemble_vacf_analysis: ensemble-averaged VACF (raw + normalised)
"""

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit
from typing import Optional, Tuple

from analysis.types import PositionsDict


# =========================================================
# Velocity (canonical forward difference)
# =========================================================

def compute_velocity(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]:
    """
    Compute forward-difference velocity from a position sequence.

    Parameters
    ----------
    r : ndarray, shape (N, 2)
        Ordered (x, y) positions.
    dt : float
        Time step between consecutive frames (seconds or simulation units).

    Returns
    -------
    v : ndarray, shape (N-1, 2)
        Velocity at each inter-frame interval:
        ``v[i] = (r[i+1] - r[i]) / dt``.

    Notes
    -----
    The returned array has one fewer row than ``r``.
    This is the sole velocity definition used throughout the analysis package.
    """
    return (r[1:] - r[:-1]) / dt


# =========================================================
# Single-droplet VACF
# =========================================================

def compute_vacf(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]:
    """
    Compute the velocity autocorrelation function for a single trajectory.

    Parameters
    ----------
    r : ndarray, shape (N, 2)
        Ordered (x, y) positions.
    dt : float
        Time step between consecutive frames.

    Returns
    -------
    vacf : ndarray, shape (N-1,)
        VACF values at lags τ = 0, 1, …, N-2 (in frame units).
        ``vacf[0] = <v·v>`` (mean kinetic energy per unit mass × 2).

    Notes
    -----
    Velocity is computed with ``compute_velocity`` (forward difference).
    The dot-product estimator is used:
    ``VACF(τ) = mean_t [vx(t)·vx(t+τ) + vy(t)·vy(t+τ)]``.
    """
    v = compute_velocity(r, dt)
    n = len(v)
    vacf = np.zeros(n)
    for tau in range(n):
        corr = v[:n - tau, 0] * v[tau:, 0] + v[:n - tau, 1] * v[tau:, 1]
        vacf[tau] = np.mean(corr)
    return vacf


# =========================================================
# Exponential fit → persistence time
# =========================================================

def fit_vacf_exponential(
    tau: NDArray[np.floating],
    vacf: NDArray[np.floating],
    max_fit_fraction: Optional[float] = None,
) -> Tuple[float, float, NDArray[np.floating], NDArray[np.floating]]:
    """
    Fit the VACF to ``C(t) = A exp(-t / τ_p)`` and return the velocity
    persistence time τ_p.

    Parameters
    ----------
    tau : ndarray
        Lag times (seconds or frame units — consistent with ``vacf``).
    vacf : ndarray
        VACF values (raw or normalised; must be positive at τ=0).
    max_fit_fraction : float or None
        If set, restricts the fit to the first ``floor(max_fit_fraction * len(vacf))``
        points *before* the first zero crossing.

    Returns
    -------
    tau_p : float
        Fitted persistence time (same units as ``tau``).
    tau_p_err : float
        Standard error of τ_p from the covariance matrix.
    tau_fit : ndarray
        Lag values used in the fit.
    fit_curve : ndarray
        Fitted curve evaluated at ``tau_fit``.

    Notes
    -----
    The fit stops at the first non-positive VACF value to avoid fitting
    oscillatory or noisy tails. Returns NaN when fewer than 5 points are
    available for fitting.
    """
    def _model(t: NDArray, A: float, tau_p: float) -> NDArray:
        return A * np.exp(-t / tau_p)

    negative = np.where(vacf <= 0)[0]
    stop = negative[0] if len(negative) > 0 else len(vacf)

    if max_fit_fraction is not None:
        stop = min(stop, max(1, int(max_fit_fraction * len(vacf))))

    tau_fit = tau[:stop]
    vacf_fit = vacf[:stop]

    if len(tau_fit) < 5:
        return np.nan, np.nan, tau_fit, np.array([])

    try:
        popt, pcov = curve_fit(
            _model,
            tau_fit,
            vacf_fit,
            p0=[vacf_fit[0], float(tau_fit.max()) / 10.0],
        )
    except RuntimeError:
        return np.nan, np.nan, tau_fit, np.array([])

    tau_p = float(popt[1])
    tau_p_err = float(np.sqrt(pcov[1, 1]))
    fit_curve = _model(tau_fit, *popt)

    return tau_p, tau_p_err, tau_fit, fit_curve


# =========================================================
# Ensemble VACF
# =========================================================

def ensemble_vacf_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray[np.int_], NDArray[np.floating], NDArray[np.floating]]:
    """
    Compute the ensemble-averaged VACF (raw and normalised).

    Parameters
    ----------
    positions_dict : PositionsDict
        Mapping droplet_id → (N, 2) position array.
    dt : float
        Time step between consecutive frames.

    Returns
    -------
    tau : ndarray, shape (L,)
        Lag times in frame units: 0, 1, …, L-1.
    vacf_raw : ndarray, shape (L,)
        Ensemble-averaged raw VACF.
    vacf_norm : ndarray, shape (L,)
        VACF normalised by its value at τ=0 (i.e. ``vacf_raw / vacf_raw[0]``).
        Equals ``vacf_raw`` if ``vacf_raw[0] == 0``.

    Notes
    -----
    All per-droplet VACFs are truncated to the shortest one before averaging.
    Velocity is computed with ``compute_velocity`` (forward difference).
    """
    all_vacf = []
    for r in positions_dict.values():
        all_vacf.append(compute_vacf(r, dt))

    min_len = min(len(v) for v in all_vacf)
    ensemble_vacf = np.mean([v[:min_len] for v in all_vacf], axis=0)
    tau = np.arange(min_len)

    vacf_norm = (
        ensemble_vacf / ensemble_vacf[0]
        if ensemble_vacf[0] != 0
        else ensemble_vacf.copy()
    )

    return tau, ensemble_vacf, vacf_norm
