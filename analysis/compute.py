"""
analysis.compute
================
All physics / statistics computations on droplet trajectories.

Contents
--------
Velocity
  - cal_velocity               : finite-difference velocity
Speed
  - ensemble_speed_analysis    : ensemble-averaged speed vs time
Orientation
  - orientation_analysis       : unwrapped heading angle from velocity
MSD
  - compute_msd                : single-droplet mean squared displacement
  - ensemble_msd_analysis      : ensemble-averaged MSD + power-law exponent
VACF
  - compute_vacf               : single-droplet velocity autocorrelation
  - ensemble_vacf_analysis     : ensemble-averaged VACF (raw + normalised)
Orientation correlation
  - compute_orientation_corr           : single-droplet cos(Δθ) correlation
  - ensemble_orientation_corr_analysis : ensemble-averaged orientation corr.
PSD
  - compute_psd_components     : Welch PSD of vx and vy
  - fit_psd_powerlaw           : power-law fit to PSD
  - ensemble_psd_analysis      : ensemble-averaged PSD
"""

import numpy as np
from numpy.typing import NDArray
from typing import Dict, Tuple
from scipy.signal import welch


# =========================================================
# VELOCITY  (simple forward difference, used by analysis)
# =========================================================

def cal_velocity(
    r: NDArray[np.floating],
    dt: float = 1.0
) -> NDArray[np.floating]:
    """Forward-difference velocity: v[i] = (r[i+1] - r[i]) / dt."""
    return (r[1:] - r[:-1]) / dt


# =========================================================
# SPEED ANALYSIS
# =========================================================

def ensemble_speed_analysis(
    positions_dict: Dict[int, NDArray[np.floating]],
    dt: float = 1.0
) -> Tuple[
    NDArray[np.int_],
    NDArray[np.floating],
    NDArray[np.floating]
]:
    """
    Compute the ensemble-averaged speed and velocity vector
    across all droplets.

    Returns (time_indices, mean_speed, mean_velocity).
    """

    ids = list(positions_dict.keys())

    min_len = min(
        len(cal_velocity(positions_dict[i], dt))
        for i in ids
    )

    ensemble_speed = np.zeros(min_len)
    ensemble_v = np.zeros((min_len, 2))

    for droplet_id in ids:

        r = positions_dict[droplet_id]
        v = cal_velocity(r, dt)[:min_len]

        speed = np.sqrt(v[:, 0]**2 + v[:, 1]**2)

        ensemble_speed += speed
        ensemble_v += v

    ensemble_speed /= len(ids)
    ensemble_v /= len(ids)

    t = np.arange(min_len)

    return t, ensemble_speed, ensemble_v


# =========================================================
# ORIENTATION
# =========================================================

def orientation_analysis(
    ensemble_v: NDArray[np.floating]
) -> NDArray[np.floating]:
    """Compute the unwrapped heading angle from the ensemble velocity."""

    theta = np.arctan2(
        ensemble_v[:, 1],
        ensemble_v[:, 0]
    )

    theta = np.unwrap(theta)

    return theta


# =========================================================
# MSD  (mean squared displacement)
# =========================================================

def compute_msd(
    r: NDArray[np.floating]
) -> NDArray[np.floating]:
    """Single-droplet MSD as a function of lag time τ."""

    n = len(r)
    msd = np.zeros(n - 1)

    for tau in range(1, n):

        disp = r[tau:] - r[:-tau]
        sq_disp = np.sum(disp**2, axis=1)
        msd[tau - 1] = np.mean(sq_disp)

    return msd


def ensemble_msd_analysis(
    positions_dict: Dict[int, NDArray[np.floating]]
) -> Tuple[
    NDArray[np.int_],
    NDArray[np.floating],
    float
]:
    """
    Ensemble-averaged MSD across all droplets with a log–log
    power-law fit to extract the diffusion exponent α.

    Returns (lag_times, ensemble_msd, alpha).
    """

    all_msds = []
    ids = list(positions_dict.keys())

    for droplet_id in ids:

        r = positions_dict[droplet_id]
        msd = compute_msd(r)
        all_msds.append(msd)

    min_len = min(len(m) for m in all_msds)
    all_msds = [m[:min_len] for m in all_msds]

    ensemble_msd = np.mean(all_msds, axis=0)

    tau = np.arange(1, min_len + 1)

    coeffs = np.polyfit(
        np.log(tau),
        np.log(ensemble_msd),
        1
    )

    alpha = coeffs[0]

    return tau, ensemble_msd, alpha


# =========================================================
# VACF  (velocity autocorrelation function)
# =========================================================

def compute_vacf(
    r: NDArray[np.floating],
    dt: float = 1.0
) -> NDArray[np.floating]:
    """Single-droplet velocity autocorrelation C_v(τ)."""

    v = cal_velocity(r, dt)
    n = len(v)
    vacf = np.zeros(n)

    for tau in range(n):

        corr = (
            v[:n - tau, 0] * v[tau:, 0]
            +
            v[:n - tau, 1] * v[tau:, 1]
        )

        vacf[tau] = np.mean(corr)

    return vacf


def ensemble_vacf_analysis(
    positions_dict: Dict[int, NDArray[np.floating]],
    dt: float = 1.0
) -> Tuple[
    NDArray[np.int_],
    NDArray[np.floating],
    NDArray[np.floating]
]:
    """
    Ensemble-averaged VACF (raw and normalised by C_v(0)).

    Returns (lag_times, raw_vacf, normalised_vacf).
    """

    all_vacf = []
    ids = list(positions_dict.keys())

    for droplet_id in ids:

        r = positions_dict[droplet_id]
        vacf = compute_vacf(r, dt)
        all_vacf.append(vacf)

    min_len = min(len(v) for v in all_vacf)
    all_vacf = [v[:min_len] for v in all_vacf]

    ensemble_vacf = np.mean(all_vacf, axis=0)

    tau = np.arange(min_len)

    if ensemble_vacf[0] != 0:
        vacf_norm = (
            ensemble_vacf / ensemble_vacf[0]
        )
    else:
        vacf_norm = ensemble_vacf

    return tau, ensemble_vacf, vacf_norm


# =========================================================
# ORIENTATION CORRELATION
# =========================================================

def compute_orientation_corr(
    theta: NDArray[np.floating]
) -> NDArray[np.floating]:
    """Single-droplet ⟨cos(Δθ(τ))⟩ correlation."""

    n = len(theta)
    corr = np.zeros(n)

    for tau in range(n):

        corr_tau = np.cos(
            theta[tau:] - theta[:n - tau]
        )

        corr[tau] = np.mean(corr_tau)

    return corr


def ensemble_orientation_corr_analysis(
    positions_dict: Dict[int, NDArray[np.floating]],
    dt: float = 1.0
) -> Tuple[
    NDArray[np.int_],
    NDArray[np.floating]
]:
    """
    Ensemble-averaged orientation correlation.

    Returns (lag_times, ensemble_correlation).
    """

    all_corr = []
    ids = list(positions_dict.keys())

    for droplet_id in ids:

        r = positions_dict[droplet_id]
        v = cal_velocity(r, dt)

        theta = np.arctan2(
            v[:, 1],
            v[:, 0]
        )

        theta = np.unwrap(theta)

        corr = compute_orientation_corr(theta)
        all_corr.append(corr)

    min_len = min(len(c) for c in all_corr)
    all_corr = [c[:min_len] for c in all_corr]

    ensemble_corr = np.mean(all_corr, axis=0)

    tau = np.arange(min_len)

    return tau, ensemble_corr


# =========================================================
# PSD  (power spectral density via Welch's method)
# =========================================================

def compute_psd_components(
    r: NDArray[np.floating],
    dt: float = 1.0
):
    """Welch PSD of the vx and vy velocity components."""

    v = cal_velocity(r, dt)

    vx = v[:, 0]
    vy = v[:, 1]

    vx = vx - np.mean(vx)
    vy = vy - np.mean(vy)

    if len(vx) < 8:
        return (
            np.array([]),
            np.array([]),
            np.array([])
        )

    fs = 1.0 / dt

    nperseg = min(256, len(vx))

    freq, psd_vx = welch(
        vx,
        fs=fs,
        nperseg=nperseg
    )

    _, psd_vy = welch(
        vy,
        fs=fs,
        nperseg=nperseg
    )

    return freq, psd_vx, psd_vy


# =========================================================
# PSD POWER LAW FIT
# =========================================================

def fit_psd_powerlaw(
    freqs: NDArray[np.floating],
    psd: NDArray[np.floating],
    fmin: float = 0.1,
    fmax: float = 2.0
) -> Tuple[float, NDArray[np.floating]]:
    """
    Fit PSD ~ f^{-β} in the band [fmin, fmax].

    Returns (beta, fit_line).
    """

    mask = (
        (freqs >= fmin)
        &
        (freqs <= fmax)
        &
        (psd > 0)
    )

    fit_freqs = freqs[mask]
    fit_psd = psd[mask]

    logf = np.log10(fit_freqs)
    logp = np.log10(fit_psd)

    coeffs = np.polyfit(
        logf,
        logp,
        1
    )

    slope = coeffs[0]
    beta = -slope

    fit_line = 10**(
        coeffs[1] + coeffs[0] * logf
    )

    return beta, fit_line


def ensemble_psd_analysis(
    positions_dict: Dict[int, NDArray[np.floating]],
    dt: float = 1.0
):
    """
    Ensemble-averaged PSD of vx and vy.

    Returns (frequencies, ensemble_psd_vx, ensemble_psd_vy).
    """

    all_psd_vx = []
    all_psd_vy = []

    freq_ref = None

    ids = list(positions_dict.keys())

    for droplet_id in ids:

        r = positions_dict[droplet_id]

        freq, psd_vx, psd_vy = (
            compute_psd_components(r, dt)
        )

        if len(freq) == 0:
            continue

        all_psd_vx.append(psd_vx)
        all_psd_vy.append(psd_vy)

        if freq_ref is None:
            freq_ref = freq

    if len(all_psd_vx) == 0:
        return (
            np.array([]),
            np.array([]),
            np.array([])
        )

    min_len = min(len(p) for p in all_psd_vx)

    all_psd_vx = [
        p[:min_len]
        for p in all_psd_vx
    ]

    all_psd_vy = [
        p[:min_len]
        for p in all_psd_vy
    ]

    freq_ref = freq_ref[:min_len]

    ensemble_psd_vx = np.mean(
        all_psd_vx,
        axis=0
    )

    ensemble_psd_vy = np.mean(
        all_psd_vy,
        axis=0
    )

    return (
        freq_ref,
        ensemble_psd_vx,
        ensemble_psd_vy
    )
