"""
analysis.psd
============
Canonical power spectral density (PSD) computations.

PSD is computed via Welch's method on the mean-subtracted velocity
components vx and vy. Frequency units are always cycles per unit of ``dt``
(Hz when ``dt`` is in seconds).

All frequency axes are returned in physical units (1/dt); callers need
not rescale.

Contents
--------
- compute_psd_components : Welch PSD of vx and vy for a single trajectory
- fit_psd_powerlaw       : power-law fit PSD ~ f^{-β}
- ensemble_psd_analysis  : ensemble-averaged PSD
"""

import numpy as np
from numpy.typing import NDArray
from scipy.signal import welch
from typing import Tuple

from analysis.types import PositionsDict
from analysis.vacf import compute_velocity

try:
    import statsmodels.api as sm
except Exception:
    sm = None


# =========================================================
# Internal OLS helper (local copy to avoid circular import)
# =========================================================

def _loglog_ols(
    log_x: NDArray[np.floating],
    log_y: NDArray[np.floating],
) -> Tuple[float, float, float, NDArray[np.floating]]:
    if len(log_x) < 2:
        return np.nan, np.nan, np.nan, np.array([])
    if sm is not None:
        X = sm.add_constant(log_x)
        model = sm.OLS(log_y, X).fit()
        slope = float(model.params[1])
        slope_err = float(model.bse[1])
        r2 = float(model.rsquared)
        return slope, slope_err, r2, 10.0 ** model.predict(X)
    x_mean = np.mean(log_x)
    y_mean = np.mean(log_y)
    ss_x = np.sum((log_x - x_mean) ** 2)
    if ss_x == 0:
        return np.nan, np.nan, np.nan, np.full_like(log_x, np.nan)
    slope = np.sum((log_x - x_mean) * (log_y - y_mean)) / ss_x
    intercept = y_mean - slope * x_mean
    fitted_log_y = intercept + slope * log_x
    residuals = log_y - fitted_log_y
    rss = np.sum(residuals ** 2)
    tss = np.sum((log_y - y_mean) ** 2)
    dof = max(len(log_x) - 2, 1)
    slope_err = np.sqrt(rss / dof / ss_x)
    r2 = 1.0 - (rss / tss) if tss > 0 else np.nan
    return slope, slope_err, r2, 10.0 ** fitted_log_y


# =========================================================
# Single-trajectory PSD
# =========================================================

def compute_psd_components(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """
    Compute the Welch power spectral density of vx and vy for one trajectory.

    Parameters
    ----------
    r : ndarray, shape (N, 2)
        Ordered (x, y) positions.
    dt : float
        Time step between consecutive frames (seconds or simulation units).

    Returns
    -------
    freq : ndarray
        Frequencies in units of 1/dt (Hz when dt is in seconds).
        Empty array if the trajectory is too short.
    psd_vx : ndarray
        Welch PSD of the mean-subtracted vx component.
    psd_vy : ndarray
        Welch PSD of the mean-subtracted vy component.

    Notes
    -----
    Velocity is computed with ``compute_velocity`` (forward difference).
    Mean is subtracted from each component before computing the PSD.
    The Welch segment length is ``min(256, len(v))``.
    Returns three empty arrays if the velocity sequence has fewer than 8 points.
    """
    v = compute_velocity(r, dt)
    vx = v[:, 0] - np.mean(v[:, 0])
    vy = v[:, 1] - np.mean(v[:, 1])

    if len(vx) < 8:
        return np.array([]), np.array([]), np.array([])

    fs = 1.0 / dt
    nperseg = min(256, len(vx))
    freq, psd_vx = welch(vx, fs=fs, nperseg=nperseg)
    _, psd_vy = welch(vy, fs=fs, nperseg=nperseg)

    return freq, psd_vx, psd_vy


# =========================================================
# Power-law fit
# =========================================================

def fit_psd_powerlaw(
    freqs: NDArray[np.floating],
    psd: NDArray[np.floating],
    fmin: float = 0.1,
    fmax: float = 2.0,
) -> Tuple[float, float, float, NDArray[np.floating], NDArray[np.floating]]:
    """
    Fit PSD ~ f^{-β} on a log-log scale over [fmin, fmax].

    Parameters
    ----------
    freqs : ndarray
        Frequency values (same units as ``compute_psd_components`` output).
    psd : ndarray
        PSD values corresponding to ``freqs``.
    fmin : float
        Lower frequency bound for the fit (same units as ``freqs``).
    fmax : float
        Upper frequency bound for the fit (same units as ``freqs``).

    Returns
    -------
    beta : float
        Power-law exponent (positive means PSD decreases with frequency).
    beta_err : float
        Standard error of β.
    r2 : float
        R² of the log-log fit.
    fit_freqs : ndarray
        Frequency values within [fmin, fmax] used for fitting.
    fit_line : ndarray
        Fitted PSD values at ``fit_freqs``.

    Notes
    -----
    The fit uses OLS on log10(f) vs log10(PSD).
    Returns NaN scalars and empty arrays when fewer than 3 points are in range.
    """
    mask = (freqs >= fmin) & (freqs <= fmax) & (psd > 0)
    fit_freqs = freqs[mask]
    fit_psd = psd[mask]

    if len(fit_freqs) < 3:
        return np.nan, np.nan, np.nan, np.array([]), np.array([])

    logf = np.log10(fit_freqs)
    logp = np.log10(fit_psd)
    slope, slope_err, r2, fit_line = _loglog_ols(logf, logp)

    return -slope, slope_err, r2, fit_freqs, fit_line


# =========================================================
# Ensemble PSD
# =========================================================

def ensemble_psd_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """
    Compute the ensemble-averaged PSD of vx and vy.

    Parameters
    ----------
    positions_dict : PositionsDict
        Mapping droplet_id → (N, 2) position array.
    dt : float
        Time step between consecutive frames.

    Returns
    -------
    freq : ndarray
        Frequencies in units of 1/dt.  Empty array if no valid PSDs.
    ensemble_psd_vx : ndarray
        Ensemble-averaged PSD of vx.
    ensemble_psd_vy : ndarray
        Ensemble-averaged PSD of vy.

    Notes
    -----
    All per-trajectory PSDs are truncated to the shortest frequency axis
    before averaging. Trajectories with fewer than 8 velocity samples are
    skipped. Velocity and PSD are computed with ``compute_psd_components``.
    """
    all_psd_vx = []
    all_psd_vy = []
    freq_ref = None

    for r in positions_dict.values():
        freq, psd_vx, psd_vy = compute_psd_components(r, dt)
        if len(freq) == 0:
            continue
        all_psd_vx.append(psd_vx)
        all_psd_vy.append(psd_vy)
        if freq_ref is None:
            freq_ref = freq

    if len(all_psd_vx) == 0:
        return np.array([]), np.array([]), np.array([])

    min_len = min(len(p) for p in all_psd_vx)
    ensemble_psd_vx = np.mean([p[:min_len] for p in all_psd_vx], axis=0)
    ensemble_psd_vy = np.mean([p[:min_len] for p in all_psd_vy], axis=0)
    freq_ref = freq_ref[:min_len]

    return freq_ref, ensemble_psd_vx, ensemble_psd_vy
