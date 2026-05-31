"""
analysis.orientation
====================
Canonical orientation correlation computations.

Orientation correlation definition (used everywhere):
    C(τ) = <cos(θ(t+τ) - θ(t))>

where θ(t) = arctan2(vy(t), vx(t)) is the heading angle, unwrapped.

All time axes use frame-index integers internally.
Convert to seconds in callers with ``tau_seconds = tau * dt``.

Contents
--------
- compute_heading_angle          : unwrapped heading angle from positions
- compute_orientation_corr       : single-droplet cos(Δθ) correlation
- fit_orientation_persistence    : exponential fit → rotational persistence time τ_r
- ensemble_orientation_corr_analysis : ensemble-averaged orientation correlation
"""

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit
from typing import Optional, Tuple

from analysis.types import PositionsDict
from analysis.vacf import compute_velocity


# =========================================================
# Heading angle
# =========================================================

def compute_heading_angle(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]:
    """
    Compute the unwrapped heading angle from a position sequence.

    Parameters
    ----------
    r : ndarray, shape (N, 2)
        Ordered (x, y) positions.
    dt : float
        Time step between consecutive frames.

    Returns
    -------
    theta : ndarray, shape (N-1,)
        Unwrapped heading angle (radians) at each inter-frame interval.
        Computed as ``arctan2(vy, vx)`` after forward-difference velocity,
        then phase-unwrapped to remove 2π discontinuities.

    Notes
    -----
    Velocity is computed with ``compute_velocity`` (forward difference).
    """
    v = compute_velocity(r, dt)
    theta = np.arctan2(v[:, 1], v[:, 0])
    return np.unwrap(theta)


# =========================================================
# Single-droplet orientation correlation
# =========================================================

def compute_orientation_corr(
    theta: NDArray[np.floating],
) -> NDArray[np.floating]:
    """
    Compute the orientation correlation function for a single trajectory.

    Parameters
    ----------
    theta : ndarray, shape (M,)
        Heading angles (radians), typically unwrapped.

    Returns
    -------
    corr : ndarray, shape (M,)
        ``C(τ) = <cos(θ(t+τ) - θ(t))>`` at lags τ = 0, 1, …, M-1
        (in frame units).

    Notes
    -----
    At τ=0 the result is identically 1 (cos(0) = 1).
    """
    n = len(theta)
    corr = np.zeros(n)
    for tau in range(n):
        corr[tau] = np.mean(np.cos(theta[tau:] - theta[:n - tau]))
    return corr


# =========================================================
# Exponential fit → rotational persistence time
# =========================================================

def fit_orientation_persistence(
    tau: NDArray[np.floating],
    corr: NDArray[np.floating],
    max_fit_fraction: Optional[float] = None,
) -> Tuple[float, float, float, NDArray[np.floating], NDArray[np.floating]]:
    """
    Fit the orientation correlation to ``C(t) = exp(-t / τ_r)`` and return
    the rotational persistence time τ_r.

    Parameters
    ----------
    tau : ndarray
        Lag times (seconds or frame units — consistent with ``corr``).
    corr : ndarray
        Orientation correlation values.
    max_fit_fraction : float or None
        If set, restricts the fit to the first ``floor(max_fit_fraction * len(corr))``
        points before the first zero crossing.

    Returns
    -------
    tau_r : float
        Fitted rotational persistence time (same units as ``tau``).
    tau_r_err : float
        Standard error of τ_r from the covariance matrix.
    r2 : float
        R-squared of the fit.
    tau_fit : ndarray
        Lag values used in the fit.
    fit_curve : ndarray
        Fitted curve evaluated at ``tau_fit``.

    Notes
    -----
    The model is a pure exponential with no amplitude prefactor (C(0)=1).
    Returns NaN when fewer than 5 points are available for fitting.
    """
    def _model(t: NDArray, tau_r: float) -> NDArray:
        return np.exp(-t / tau_r)

    negative = np.where(corr <= 0)[0]
    stop = negative[0] if len(negative) > 0 else len(corr)

    if max_fit_fraction is not None:
        stop = min(stop, max(1, int(max_fit_fraction * len(corr))))

    tau_fit = tau[:stop]
    corr_fit = corr[:stop]

    if len(tau_fit) < 5:
        return np.nan, np.nan, np.nan, tau_fit, np.array([])

    try:
        popt, pcov = curve_fit(
            _model,
            tau_fit,
            corr_fit,
            p0=[float(tau_fit.max()) / 10.0],
        )
    except RuntimeError:
        return np.nan, np.nan, np.nan, tau_fit, np.array([])

    tau_r = float(popt[0])
    tau_r_err = float(np.sqrt(pcov[0, 0]))
    fit_curve = _model(tau_fit, tau_r)
    
    ss_res = np.sum((corr_fit - fit_curve) ** 2)
    ss_tot = np.sum((corr_fit - np.mean(corr_fit)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    return tau_r, tau_r_err, r2, tau_fit, fit_curve


# =========================================================
# Ensemble orientation correlation
# =========================================================

def ensemble_orientation_corr_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]:
    """
    Compute the ensemble-averaged orientation correlation.

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
    ensemble_corr : ndarray, shape (L,)
        Ensemble-averaged ``<cos(θ(t+τ) - θ(t))>``.

    Notes
    -----
    Velocity and heading angle are computed with ``compute_velocity`` and
    ``compute_heading_angle`` respectively (forward difference, unwrapped).
    All per-droplet correlations are truncated to the shortest before averaging.
    """
    all_corr = []
    for r in positions_dict.values():
        theta = compute_heading_angle(r, dt)
        all_corr.append(compute_orientation_corr(theta))

    min_len = min(len(c) for c in all_corr)
    ensemble_corr = np.mean([c[:min_len] for c in all_corr], axis=0)
    tau = np.arange(min_len)

    return tau, ensemble_corr
