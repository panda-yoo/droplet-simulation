"""
analysis.msd
============
Canonical mean-squared displacement (MSD) computations.

All time axes use frame-index integers internally.
Convert to seconds in callers with ``tau_seconds = tau * dt``.

Contents
--------
- compute_msd               : single-droplet MSD vs lag (frame units)
- fit_msd_power_law         : log-log power-law fit MSD ~ τ^α
- ensemble_msd_analysis     : ensemble-averaged MSD + exponent
- compute_local_msd_exponent: d log MSD / d log τ at every lag
- compute_pair_weighted_msd : pair-weighted MSD across all tracks
"""

import numpy as np
from numpy.typing import NDArray
from typing import Dict, Optional, Tuple

try:
    import statsmodels.api as sm
except Exception:
    sm = None
    print(Exception)

from analysis.types import PositionsDict


# =========================================================
# Internal OLS helper
# =========================================================

def _loglog_ols(
    log_x: NDArray[np.floating],
    log_y: NDArray[np.floating],
) -> Tuple[float, float, float, NDArray[np.floating]]:
    """
    Ordinary-least-squares fit on log-log data.

    Parameters
    ----------
    log_x : ndarray
        log10 of the independent variable.
    log_y : ndarray
        log10 of the dependent variable.

    Returns
    -------
    slope : float
    slope_err : float
    r2 : float
    fitted_y : ndarray
        Fitted values in original (non-log) units: 10^(intercept + slope*log_x).
    """
    if len(log_x) < 2:
        return np.nan, np.nan, np.nan, np.array([])

    if sm is not None:
        X = sm.add_constant(log_x)
        model = sm.OLS(log_y, X).fit()
        slope = float(model.params[1])
        slope_err = float(model.bse[1])
        r2 = float(model.rsquared)
        fitted_y = 10.0 ** model.predict(X)
        return slope, slope_err, r2, fitted_y

    # # Fallback: manual OLS
    # x_mean = np.mean(log_x)
    # y_mean = np.mean(log_y)
    # ss_x = np.sum((log_x - x_mean) ** 2)
    # if ss_x == 0:
    #     return np.nan, np.nan, np.nan, np.full_like(log_x, np.nan)
    # slope = np.sum((log_x - x_mean) * (log_y - y_mean)) / ss_x
    # intercept = y_mean - slope * x_mean
    # fitted_log_y = intercept + slope * log_x
    # residuals = log_y - fitted_log_y
    # rss = np.sum(residuals ** 2)
    # tss = np.sum((log_y - y_mean) ** 2)
    # dof = max(len(log_x) - 2, 1)
    # sigma2 = rss / dof
    # slope_err = np.sqrt(sigma2 / ss_x)
    # r2 = 1.0 - (rss / tss) if tss > 0 else np.nan
    # return slope, slope_err, r2, 10.0 ** fitted_log_y


# =========================================================
# Single-droplet MSD
# =========================================================

def compute_msd(
    r: NDArray[np.floating],
    lag_fraction: float = 1.0,
    minimum_pairs: Optional[int] = None,
) -> NDArray[np.floating]:
    """
    Compute the mean-squared displacement for a single trajectory.

    Parameters
    ----------
    r : ndarray, shape (N, 2)
        Sequence of (x, y) positions in frame order.
    lag_fraction : float
        Fraction of ``N - 1`` to use as the maximum lag.
        Must be in ``(0, 1]``.
    minimum_pairs : int or None
        If set, the maximum lag is further capped so that at least
        ``minimum_pairs`` displacement pairs contribute to every lag bin.

    Returns
    -------
    msd : ndarray, shape (max_lag,)
        MSD values for lags τ = 1, 2, …, max_lag (in frame units).

    Notes
    -----
    Uses the direct (non-FFT) sliding-window estimator:
    ``MSD(τ) = mean_t ||r(t+τ) − r(t)||²``.
    """
    n = len(r)
    max_lag = max(1, int(lag_fraction * (n - 1)))
    if minimum_pairs is not None:
        max_lag = min(max_lag, max(1, n - minimum_pairs))

    msd = np.zeros(max_lag)
    for tau in range(1, max_lag + 1):
        disp = r[tau:] - r[:-tau]
        msd[tau - 1] = np.mean(np.sum(disp ** 2, axis=1))

    return msd


# =========================================================
# Power-law fit
# =========================================================

def fit_msd_power_law(
    msd: NDArray[np.floating],
    fit_fraction: float = 0.3,
    minimum_pairs: Optional[int] = None,
) -> Tuple[float, float, float, NDArray[np.floating]]:
    """
    Fit MSD ~ τ^α on a log-log scale.

    Parameters
    ----------
    msd : ndarray, shape (M,)
        MSD values at lags τ = 1, …, M (in frame units).
    fit_fraction : float
        Fraction of ``M`` used for fitting; the fit is restricted to
        τ ∈ [1, floor(fit_fraction * M)].
    minimum_pairs : int or None
        If set, further caps the fit range so at least ``minimum_pairs``
        displacement pairs exist at every included lag.

    Returns
    -------
    alpha : float
        Power-law exponent.
    alpha_err : float
        Standard error of alpha.
    r2 : float
        Coefficient of determination (R²) of the log-log fit.
    fit_line_full : ndarray, shape (M,)
        Fitted MSD values over the fit region; NaN outside.

    Notes
    -----
    The fit is performed on log10(τ) vs log10(MSD) using OLS.
    ``statsmodels`` is used when available; otherwise a pure-NumPy
    fallback is used (identical formulas).
    """
    m = len(msd)
    tau = np.arange(1, m + 1)

    max_fit = max(1, int(fit_fraction * m))
    if minimum_pairs is not None:
        max_fit = min(max_fit, max(1, m - minimum_pairs))

    tau_fit = tau[:max_fit]
    msd_fit = msd[:max_fit]
    mask = (tau_fit > 0) & (msd_fit > 0)

    fit_line_full = np.full(m, np.nan, dtype=float)

    if mask.sum() < 2:
        return np.nan, np.nan, np.nan, fit_line_full

    log_tau = np.log10(tau_fit[mask])
    log_msd = np.log10(msd_fit[mask])

    slope, slope_err, r2, fitted_y = _loglog_ols(log_tau, log_msd)
    fit_line_full[:len(fitted_y)] = fitted_y

    return slope, slope_err, r2, fit_line_full


# =========================================================
# Ensemble MSD + exponent
# =========================================================

def ensemble_msd_analysis(
    positions_dict: PositionsDict,
    lag_fraction: float = 0.3,
    fit_fraction: float = 0.3,
    minimum_pairs: Optional[int] = None,
) -> Tuple[
    NDArray[np.int_],
    NDArray[np.floating],
    float,
    float,
    float,
    NDArray[np.floating],
]:
    """
    Compute ensemble-averaged MSD and its power-law exponent.

    Parameters
    ----------
    positions_dict : PositionsDict
        Mapping droplet_id → (N, 2) position array.
    lag_fraction : float
        Fraction of each trajectory length used as the maximum lag.
    fit_fraction : float
        Fraction of the MSD curve included in the log-log power-law fit.
    minimum_pairs : int or None
        Minimum number of displacement pairs required at every lag.

    Returns
    -------
    tau : ndarray, shape (L,)
        Lag times in frame units: 1, 2, …, L.
    ensemble_msd : ndarray, shape (L,)
        Ensemble-averaged MSD: mean over all per-droplet MSDs.
    alpha : float
        Power-law exponent from log-log fit.
    alpha_err : float
        Standard error of alpha.
    r2 : float
        R² of the log-log fit.
    fit_line_full : ndarray, shape (L,)
        Fitted MSD values over the fit region; NaN outside.

    Notes
    -----
    All per-droplet MSDs are truncated to the length of the shortest one
    before averaging. ``fit_msd_power_law`` is called on the ensemble MSD
    using the same ``fit_fraction`` and ``minimum_pairs``; batch and
    plotting code that call this function will therefore always produce
    identical exponents.
    """
    all_msds = []
    for r in positions_dict.values():
        all_msds.append(
            compute_msd(r, lag_fraction=lag_fraction, minimum_pairs=minimum_pairs)
        )

    min_len = min(len(m) for m in all_msds)
    ensemble_msd = np.mean([m[:min_len] for m in all_msds], axis=0)
    tau = np.arange(1, min_len + 1)

    alpha, alpha_err, r2, fit_line_full = fit_msd_power_law(
        ensemble_msd,
        fit_fraction=fit_fraction,
        minimum_pairs=minimum_pairs,
    )

    return tau, ensemble_msd, alpha, alpha_err, r2, fit_line_full


# =========================================================
# Local MSD exponent
# =========================================================

def compute_local_msd_exponent(
    tau: NDArray[np.floating],
    msd: NDArray[np.floating],
) -> NDArray[np.floating]:
    """
    Compute the local (scale-dependent) MSD exponent α(τ).

    Parameters
    ----------
    tau : ndarray, shape (L,)
        Lag times (any consistent units; must be strictly positive and
        strictly increasing).
    msd : ndarray, shape (L,)
        MSD values corresponding to each lag (must be strictly positive).

    Returns
    -------
    alpha_local : ndarray, shape (L,)
        Local exponent evaluated as d log10(MSD) / d log10(τ) via
        ``numpy.gradient`` (central differences interior, one-sided at edges).

    Notes
    -----
    ``tau`` may be in frame units or seconds; the gradient is
    unit-agnostic because both axes are log10-transformed first.
    Raise ``ValueError`` if inputs are inconsistent.
    """
    tau = np.asarray(tau, dtype=float)
    msd = np.asarray(msd, dtype=float)

    if tau.ndim != 1 or msd.ndim != 1:
        raise ValueError("tau and msd must be 1-D arrays.")
    if len(tau) != len(msd):
        raise ValueError("tau and msd must have the same length.")
    if len(tau) < 3:
        raise ValueError("Need at least 3 points to compute a gradient.")
    if not np.all(np.isfinite(tau)) or not np.all(np.isfinite(msd)):
        raise ValueError("tau and msd must contain only finite values.")
    if np.any(tau <= 0) or np.any(msd <= 0):
        raise ValueError("tau and msd must be strictly positive.")
    if not np.all(np.diff(tau) > 0):
        raise ValueError("tau must be strictly increasing.")

    return np.gradient(np.log10(msd), np.log10(tau))


# =========================================================
# Pair-weighted MSD
# =========================================================

def compute_pair_weighted_msd(
    positions_dict: PositionsDict,
    lag_fraction: float = 0.3,
    minimum_pairs: Optional[int] = None,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]:
    """
    Compute pair-weighted MSD across all trajectories (including partial tracks).

    Parameters
    ----------
    positions_dict : PositionsDict
        Mapping droplet_id → (N_i, 2) position array.
        Trajectories need not have equal length.
    lag_fraction : float
        Fraction of the shortest trajectory used as the maximum lag.
    minimum_pairs : int or None
        Minimum number of displacement pairs required at every lag.

    Returns
    -------
    tau : ndarray, shape (L,)
        Lag times in frame units: 1, 2, …, L.
    msd_weighted : ndarray, shape (L,)
        Pair-weighted MSD: total squared displacement / total pair count at each lag.

    Notes
    -----
    Unlike ``ensemble_msd_analysis``, which averages per-trajectory MSDs,
    this estimator weights each displacement pair equally regardless of
    which trajectory it came from. It is therefore suitable for diagnosing
    selection bias when trajectories have unequal lengths.
    """
    if len(positions_dict) == 0:
        return np.array([], dtype=int), np.array([], dtype=float)

    lengths = [len(r) for r in positions_dict.values()]
    min_len = min(lengths)
    max_lag = max(1, int(lag_fraction * (min_len - 1)))
    if minimum_pairs is not None:
        max_lag = min(max_lag, max(1, min_len - minimum_pairs))

    msd_sum = np.zeros(max_lag, dtype=float)
    pair_counts = np.zeros(max_lag, dtype=float)

    for r in positions_dict.values():
        n = len(r)
        for tau in range(1, max_lag + 1):
            if n <= tau:
                break
            disp = r[tau:] - r[:-tau]
            sq = np.sum(disp ** 2, axis=1)
            msd_sum[tau - 1] += np.sum(sq)
            pair_counts[tau - 1] += len(sq)

    with np.errstate(invalid="ignore", divide="ignore"):
        msd_weighted = msd_sum / pair_counts

    return np.arange(1, max_lag + 1), msd_weighted
