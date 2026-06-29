"""
analysis.plotting
=================
All plotting routines for the droplet analysis pipeline.

Every plot function accepts pre-computed quantities and a label string.
No physics is computed here; all quantities are passed in from callers
that use the canonical functions in msd.py, vacf.py, orientation.py, psd.py.

Time convention
---------------
- Internal computations use frame-index lags (integers).
- Every plotting function receives ``tau_seconds`` (= tau * dt) for axis
  labels and ``fit_line`` arrays already aligned to ``tau_seconds``.

Contents
--------
- plot_trajectories
- plot_msd_loglog
- plot_vacf
- plot_orientation_corr
- plot_psd_component
- plot_local_msd_exponent
- plot_msd_comparison
- plot_track_length_hist
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from analysis.plot_style import set_publication_style, StyleTokens, apply_grid, get_droplet_color, apply_legend, add_diffusive_line, add_ballistic_line, add_zero_line

# Apply global publication style
set_publication_style()
from numpy.typing import NDArray
from typing import Optional

from analysis.types import PositionsDict
from seaborn import heatmap

# =========================================================
# Trajectories
# =========================================================

def plot_trajectories(
    data: PositionsDict,
    label: str,
    output_dir: str,
    total_paths: int,
) -> None:
    """
    Plot consistent (full-length) droplet trajectories in the x-y plane.

    Parameters
    ----------
    data : PositionsDict
        Mapping droplet_id → (N, 2) position array (full-length tracks only).
    label : str
        Dataset label used in the title and output filename.
    output_dir : str
        Directory where the figure is saved.
    total_paths : int
        Total number of droplets detected (including discarded tracks).
    """
    n_ensemble = len(data)

    fig, ax = plt.subplots(figsize=(7, 7))

    for idx, (droplet_id, pos) in enumerate(data.items()):
        c = get_droplet_color(droplet_id)
        ax.plot(pos[:, 0], pos[:, 1], linewidth=1.5, color=c,
                label=f"Droplet {droplet_id}")
        ax.plot(pos[0, 0], pos[0, 1], "o", color=c, markersize=6)
        ax.plot(pos[-1, 0], pos[-1, 1], "s", color=c, markersize=6)

    ax.set_title(
        f"{label} : Trajectories\n"
        f"N = {n_ensemble} droplets"
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    apply_legend(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{label}_trajectories.png"), dpi=300)
    plt.close(fig)


# =========================================================
# Log-log MSD
# =========================================================

def plot_msd_loglog(
    tau_seconds: NDArray[np.floating],
    ensemble_msd: NDArray[np.floating],
    alpha: float,
    alpha_err: float,
    r2_msd: float,
    fit_line_full: NDArray[np.floating],
    label: str,
    output_dir: str,
    n_ensemble: int,
) -> None:
    """
    Plot ensemble MSD on log-log axes with the power-law fit overlay.

    Parameters
    ----------
    tau_seconds : ndarray
        Lag times in seconds (tau * dt).
    ensemble_msd : ndarray
        Ensemble-averaged MSD.
    alpha : float
        Power-law exponent from ``fit_msd_power_law``.
    fit_line_full : ndarray
        Fitted MSD values (NaN outside fit region), from ``fit_msd_power_law``.
    label : str
        Dataset label.
    output_dir : str
        Directory where the figure is saved.
    n_ensemble : int
        Number of droplets in the ensemble.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(tau_seconds, ensemble_msd, color=StyleTokens.MSD, label="Ensemble MSD")

    fit_mask = ~np.isnan(fit_line_full)
    if np.any(fit_mask):
        ax.loglog(
            tau_seconds[fit_mask],
            fit_line_full[fit_mask],
            "k:",
            linewidth=2.5,
            label=rf"Fit ($\alpha = {alpha:.2f} \pm {alpha_err:.2f}$, $R^2 = {r2_msd:.2f}$)",
        )

    ax.set_title(
        f"{label} : MSD\n"
        rf"$\alpha = {alpha:.2f} \pm {alpha_err:.2f}$ | N = {n_ensemble} droplets"
    )
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"$\langle \mathrm{MSD}(\tau)\rangle$")
    apply_legend(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{label}_loglog_msd.png"), dpi=300)
    plt.close(fig)


# =========================================================
# VACF
# =========================================================

def plot_vacf(
    tau_vacf_seconds: NDArray[np.floating],
    vacf_norm: NDArray[np.floating],
    label: str,
    output_dir: str,
    n_ensemble: int,
    tau_p_fit: Optional[NDArray[np.floating]] = None,
    vacf_fit_curve: Optional[NDArray[np.floating]] = None,
    tau_p: Optional[float] = None,
    tau_p_err: Optional[float] = None,
    r2: float = 0.0,
) -> None:
    """
    Plot normalised VACF (first quarter of the lag range).

    Parameters
    ----------
    tau_vacf_seconds : ndarray
        Lag times in seconds.
    vacf_norm : ndarray
        Normalised VACF (divided by its value at τ=0).
    label : str
        Dataset label.
    output_dir : str
        Directory where the figure is saved.
    n_ensemble : int
        Number of droplets in the ensemble.
    tau_p_fit : ndarray or None
        Lag values used in the exponential fit (seconds).
    vacf_fit_curve : ndarray or None
        Fitted curve values at ``tau_p_fit``.
    tau_p : float or None
        Fitted persistence time (seconds); shown in legend when not NaN.
    """
    cut = len(tau_vacf_seconds) // 4 or len(tau_vacf_seconds)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(tau_vacf_seconds[:cut], vacf_norm[:cut], color=StyleTokens.VACF, label=r"$C_v(\tau)/C_v(0)$")

    if (
        tau_p_fit is not None
        and vacf_fit_curve is not None
        and len(vacf_fit_curve) > 0
        and tau_p is not None
        and not np.isnan(tau_p)
    ):
        mask = tau_p_fit <= tau_vacf_seconds[cut - 1]
        if np.any(mask):
            err_str = rf" \pm {tau_p_err:.1f}" if tau_p_err is not None else ""
            legend_label = rf"Fit ($\tau_p = {tau_p:.1f}{err_str}\,\mathrm{{s}}$, $R^2 = {r2:.2f}$)"
            ax.plot(tau_p_fit[mask], vacf_fit_curve[mask], "k--",
                    linewidth=2.5, label=legend_label)

    ax.set_title(
        f"{label} : VACF\n"
        f"N = {n_ensemble} droplets"
    )
    add_zero_line(ax)
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"$C_v(\tau)/C_v(0)$")
    apply_legend(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{label}_vacf.png"), dpi=300)
    plt.close(fig)


# =========================================================
# Orientation correlation
# =========================================================

def plot_orientation_corr(
    tau_orient_seconds: NDArray[np.floating],
    orient_corr: NDArray[np.floating],
    label: str,
    output_dir: str,
    n_ensemble: int,
    tau_r_fit: Optional[NDArray[np.floating]] = None,
    orient_fit_curve: Optional[NDArray[np.floating]] = None,
    tau_r: Optional[float] = None,
    tau_r_err: Optional[float] = None,
    r2: float = 0.0,
) -> None:
    """
    Plot ensemble orientation correlation.

    Parameters
    ----------
    tau_orient_seconds : ndarray
        Lag times in seconds.
    orient_corr : ndarray
        Ensemble-averaged ``<cos(Δθ)>``.
    label : str
        Dataset label.
    output_dir : str
        Directory where the figure is saved.
    n_ensemble : int
        Number of droplets in the ensemble.
    tau_r_fit : ndarray or None
        Lag values used in the exponential fit (seconds).
    orient_fit_curve : ndarray or None
        Fitted curve values at ``tau_r_fit``.
    tau_r : float or None
        Fitted rotational persistence time (seconds); shown in legend when not NaN.
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(tau_orient_seconds, orient_corr, color=StyleTokens.ORIENT_CORR,
            label=r"$\langle\cos(\Delta\theta)\rangle$")

    if (
        tau_r_fit is not None
        and orient_fit_curve is not None
        and len(orient_fit_curve) > 0
        and tau_r is not None
        and not np.isnan(tau_r)
    ):
        err_str = rf" \pm {tau_r_err:.1f}" if tau_r_err is not None else ""
        ax.plot(tau_r_fit, orient_fit_curve, "k--", linewidth=2.5,
                label=rf"Fit ($\tau_r = {tau_r:.1f}{err_str}\,\mathrm{{s}}$, $R^2 = {r2:.2f}$)")

    ax.set_title(
        f"{label} : Orientation Correlation\n"
        f"N = {n_ensemble} droplets"
    )
    add_zero_line(ax)
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"$\langle \cos(\Delta \theta)\rangle$")
    apply_legend(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{label}_orientation_corr.png"), dpi=300)
    plt.close(fig)


# =========================================================
# PSD (single component)
# =========================================================

def plot_psd_component(
    freq_psd: NDArray[np.floating],
    psd: NDArray[np.floating],
    beta: float,
    beta_err: float,
    r2: float,
    fit_freqs: NDArray[np.floating],
    fit_line: NDArray[np.floating],
    component: str,
    label: str,
    output_dir: str,
    n_ensemble: int,
) -> None:
    """
    Plot PSD of one velocity component (vx or vy) on log-log axes.

    Parameters
    ----------
    freq_psd : ndarray
        Frequency values (Hz).
    psd : ndarray
        PSD values.
    beta : float
        Power-law exponent from ``fit_psd_powerlaw``.
    fit_freqs : ndarray
        Frequency values used in the fit.
    fit_line : ndarray
        Fitted PSD values at ``fit_freqs``.
    component : str
        Either ``"vx"`` or ``"vy"``.
    label : str
        Dataset label.
    output_dir : str
        Directory where the figure is saved.
    n_ensemble : int
        Number of droplets in the ensemble.
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.loglog(freq_psd, psd, color=StyleTokens.PSD_X if component == "vx" else StyleTokens.PSD_Y, label="PSD")

    if len(fit_freqs) > 0 and not np.isnan(beta):
        ax.loglog(fit_freqs, fit_line, "k:", linewidth=2.5,
                  label=rf"Fit ($\beta = {beta:.2f} \pm {beta_err:.2f}$, $R^2 = {r2:.2f}$)")

    subscript = "x" if component == "vx" else "y"
    ax.set_title(
        f"{label} : PSD(v{subscript})\n"
        rf"N = {n_ensemble} droplets"
    )
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel(rf"$S_{{v_{subscript}}}(f)$")
    apply_legend(ax)
    fig.tight_layout()
    suffix = "vx" if component == "vx" else "vy"
    fig.savefig(os.path.join(output_dir, f"{label}_psd_{suffix}.png"), dpi=300)
    plt.close(fig)


# =========================================================
# Local MSD exponent
# =========================================================

def plot_local_msd_exponent(
    tau_seconds: NDArray[np.floating],
    d_msd: NDArray[np.floating],
    alpha: float,
    alpha_err: float,
    label: str,
    output_dir: str,
    n_ensemble: int,
    window_length: Optional[int] = None,
    polyorder: int = 2,
) -> None:
    """
    Plot the local MSD exponent α(τ) vs lag time on a semi-log x-axis.

    Parameters
    ----------
    tau_seconds : ndarray
        Lag times in seconds.
    d_msd : ndarray
        Local MSD exponent from ``compute_local_msd_exponent``.
    alpha : float
        Global power-law exponent (shown as a horizontal reference line).
    label : str
        Dataset label.
    output_dir : str
        Directory where the figure is saved.
    n_ensemble : int
        Number of droplets in the ensemble.
    window_length : int or None
        If set, Savitzky-Golay smoothing is applied with this window length.
    polyorder : int
        Polynomial order for Savitzky-Golay smoothing (used only when
        ``window_length`` is not None).
    """
    valid = slice(1, -1)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogx(tau_seconds[valid], d_msd[valid], alpha=0.6, color=StyleTokens.LOCAL_ALPHA,
                label=r"$\alpha(\tau)$ (raw)")

    if window_length is not None:
        try:
            from scipy.signal import savgol_filter
            d_smooth = savgol_filter(d_msd, window_length=window_length,
                                     polyorder=polyorder)
            ax.semilogx(tau_seconds[valid], d_smooth[valid],
                        linewidth=2.0, label=r"$\alpha(\tau)$ (smoothed)")
        except Exception:
            pass

    ax.axhline(alpha, linestyle="--", color="k", lw=2.5, label=rf"Global $\alpha = {alpha:.2f} \pm {alpha_err:.2f}$")
    add_diffusive_line(ax)
    add_ballistic_line(ax)

    ax.set_title(
        f"{label} : Local MSD Exponent\n"
        f"N = {n_ensemble} droplets"
    )
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"$\alpha(\tau)$")
    apply_legend(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{label}_alpha_vs_tau.png"), dpi=300)
    plt.close(fig)


# =========================================================
# MSD comparison (ensemble vs pair-weighted)
# =========================================================

def plot_msd_comparison(
    tau_seconds: NDArray[np.floating],
    ensemble_msd: NDArray[np.floating],
    tau_weighted_seconds: NDArray[np.floating],
    msd_weighted: NDArray[np.floating],
    label: str,
    output_dir: str,
    n_ensemble: int,
) -> None:
    """
    Plot full-length ensemble MSD vs pair-weighted MSD on log-log axes.

    Parameters
    ----------
    tau_seconds : ndarray
        Lag times for the ensemble MSD (seconds).
    ensemble_msd : ndarray
        Ensemble-averaged MSD (full-length tracks only).
    tau_weighted_seconds : ndarray
        Lag times for the pair-weighted MSD (seconds).
    msd_weighted : ndarray
        Pair-weighted MSD (all tracks).
    label : str
        Dataset label.
    output_dir : str
        Directory where the figure is saved.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(tau_seconds, ensemble_msd, color=StyleTokens.MSD, label="Full-length ensemble")
    ax.loglog(tau_weighted_seconds, msd_weighted, color=StyleTokens.ORIENT_CORR,
              label="Pair-weighted (all tracks)")
    import logging
    if len(msd_weighted) == len(ensemble_msd):
        logging.warning("Pair count matches ensemble length! Is pair_count being reported as droplets?")

    ax.set_title(f"{label} : MSD Comparison\nN = {n_ensemble} droplets")
    ax.text(0.05, 0.95, f"{len(msd_weighted)} displacement pairs", 
            transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"$\langle \mathrm{MSD}(\tau)\rangle$")
    apply_legend(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{label}_msd_pair_weighted.png"), dpi=300)
    plt.close(fig)


# =========================================================
# Track length histogram
# =========================================================

def plot_track_length_hist(
    retained_lengths: list,
    discarded_lengths: list,
    label: str,
    output_dir: str,
    n_ensemble: int,
) -> None:
    """
    Plot a histogram of retained vs discarded track lengths.

    Parameters
    ----------
    retained_lengths : list of int
        Lengths of full-length tracks.
    discarded_lengths : list of int
        Lengths of partial (discarded) tracks.
    label : str
        Dataset label.
    output_dir : str
        Directory where the figure is saved.
    """
    total = len(retained_lengths) + len(discarded_lengths)
    bins = min(20, max(5, total))

    fig, ax = plt.subplots(figsize=(6, 4))
    if retained_lengths:
        ax.hist(retained_lengths, bins=bins, color=StyleTokens.MSD, alpha=0.7, label="retained")
    if discarded_lengths:
        ax.hist(discarded_lengths, bins=bins, color=StyleTokens.END, alpha=0.7, label="discarded")
    ax.set_title(f"{label} : Track Length Distribution")
    ax.set_xlabel("Track length (frames)")
    ax.set_ylabel("Count")
    apply_legend(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{label}_track_length_hist.png"), dpi=300)
    plt.close(fig)


def plot_vacf_heatmap_btw_droplets(
    corr_matrix : NDArray[np.floating],
    droplet_ids : list,save_plot : bool = False,
    output_dir : str = "",label : str = "",
    lag_time : int = 0) -> None:
    
    fig, ax = plt.subplots(figsize=(6, 4))
  
    ax.set_xlabel(r"droplets")
    ax.set_ylabel(r"droplets")
    
    heatmap(corr_matrix[lag_time],xticklabels=droplet_ids,
            yticklabels=droplet_ids,annot=True,fmt=".2f",
            cmap="coolwarm",ax=ax)
    plt.show()
    plt.close(fig)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{label}_msd_pair_weighted.png"), dpi=300)


