"""
single_droplet_analysis.py
===========================
Publication-quality characterisation of a single droplet trajectory.

All physical quantities are computed using the canonical implementations
in the ``analysis/`` package — the same functions used by the ensemble
pipeline.  Nothing is recomputed or reimplemented here.

Typical usage
-------------
From a script::

    from single_droplet_analysis import analyze_single_droplet
    from analysis.trajectories import load_positions
    from pathlib import Path

    positions, _, _ = load_positions(Path("data.csv"))
    r = positions[0]          # (N, 2) array for one droplet

    result = analyze_single_droplet(
        positions=r,
        dt=0.2,
        droplet_id=0,
        experiment_name="for_2.1_droplets",
        output_dir=Path("output/plots/single_droplet"),
    )
    print(result)

From a notebook::

    from single_droplet_analysis import analyze_single_droplet
    import numpy as np
    from pathlib import Path

    result = analyze_single_droplet(
        positions=positions_dict[42],
        dt=1.0 / 30.0,
        droplet_id=42,
        experiment_name="experiment_A",
        output_dir=Path("output/plots/single_droplet"),
    )

Outputs
-------
- ``{output_dir}/{experiment_name}_droplet_{id}_summary.png``
  — 2 × 3 summary figure
- ``SingleDropletAnalysisResult`` dataclass with all scalar results
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from numpy.typing import NDArray

# ── Ensure the project root is on sys.path so that imports work whether
#    this module is run directly (python single_droplet_analysis.py) or
#    imported from a notebook sitting in a sub-directory.
_PROJECT_ROOT = Path(__file__).resolve().parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# ── Canonical project imports ────────────────────────────────────────────────
from analysis.vacf import compute_velocity, compute_vacf, fit_vacf_exponential
from analysis.msd import compute_msd, fit_msd_power_law, compute_local_msd_exponent
from analysis.orientation import (
    compute_heading_angle,
    compute_orientation_corr,
    fit_orientation_persistence,
)
from analysis.psd import compute_psd_components, fit_psd_powerlaw
from analysis.trajectories import load_positions

# ── Typography ───────────────────────────────────────────────────────────────
from analysis.plot_style import set_publication_style, StyleTokens, apply_grid, get_droplet_color
set_publication_style()

__all__ = ["SingleDropletAnalysisResult", "analyze_single_droplet", "analyze_all_droplets"]


# ============================================================================
# Result dataclass
# ============================================================================


@dataclass
class SingleDropletAnalysisResult:
    """
    Scalar summary of a single-droplet analysis.

    All time-valued fields are in seconds (``dt``-scaled).
    Speed fields are in the same position units per second as the input.

    Attributes
    ----------
    droplet_id : int or str
        Identifier of the analysed droplet.
    experiment_name : str
        Label of the experiment or dataset.
    track_length : int
        Number of frames N in the trajectory.
    mean_speed : float
        Time-averaged speed ``⟨|v|⟩`` (position units / s).
    rms_speed : float
        Root-mean-squared speed ``√⟨|v|²⟩`` (position units / s).
    msd_alpha : float
        Power-law exponent α from ``MSD(τ) ~ τ^α``.
    msd_alpha_error : float
        Standard error of α (OLS).
    msd_r2 : float
        R² of the log-log power-law fit.
    persistence_time : float
        Velocity persistence time τ_p (seconds); NaN if fit failed.
    rotational_persistence_time : float
        Rotational persistence time τ_r (seconds); NaN if fit failed.
    psd_beta_x : float
        PSD power-law exponent β for vx; NaN if fit failed.
    psd_beta_y : float
        PSD power-law exponent β for vy; NaN if fit failed.
    turning_angle_mean : float
        Mean turning angle Δθ (radians).
    turning_angle_std : float
        Standard deviation of turning angles (radians).
    figure_path : Path or None
        Absolute path to the saved summary figure.
    """

    droplet_id: int | str
    experiment_name: str
    track_length: int
    mean_speed: float
    rms_speed: float
    msd_alpha: float
    msd_alpha_error: float
    msd_r2: float
    persistence_time: float
    rotational_persistence_time: float
    psd_beta_x: float
    psd_beta_y: float
    turning_angle_mean: float
    turning_angle_std: float
    figure_path: Optional[Path] = field(default=None)

    def __str__(self) -> str:
        lines = [
            f"=== SingleDropletAnalysisResult ===",
            f"  Experiment     : {self.experiment_name}",
            f"  Droplet ID     : {self.droplet_id}",
            f"  Track length   : {self.track_length} frames",
            f"",
            f"  Speed          : mean = {self.mean_speed:.4g}  rms = {self.rms_speed:.4g}",
            f"  MSD exponent α : {self.msd_alpha:.4f} ± {self.msd_alpha_error:.4f}  R² = {self.msd_r2:.4f}",
            f"  Persist. τ_p   : {self.persistence_time:.4g} s",
            f"  Rotat.   τ_r   : {self.rotational_persistence_time:.4g} s",
            f"  PSD β_x        : {self.psd_beta_x:.4f}",
            f"  PSD β_y        : {self.psd_beta_y:.4f}",
            f"  Turn. Δθ       : mean = {np.degrees(self.turning_angle_mean):.2f}°"
            f"  std = {np.degrees(self.turning_angle_std):.2f}°",
            f"",
            f"  Figure         : {self.figure_path}",
        ]
        return "\n".join(lines)


# ============================================================================
# Internal helpers
# ============================================================================


def _turning_angles(theta: NDArray[np.floating]) -> NDArray[np.floating]:
    """
    Compute successive turning angles from an unwrapped heading sequence.

    Parameters
    ----------
    theta : ndarray, shape (M,)
        Unwrapped heading angles in radians (from ``compute_heading_angle``).

    Returns
    -------
    delta_theta : ndarray, shape (M-1,)
        Successive differences ``theta[i+1] - theta[i]``, wrapped to
        ``(-π, π]``.
    """
    raw = np.diff(theta)
    # Wrap to (-pi, pi]
    return (raw + np.pi) % (2.0 * np.pi) - np.pi


def _safe_label(x: float, fmt: str = ".3f") -> str:
    """Format a scalar for a legend label; return 'N/A' on NaN."""
    return f"{x:{fmt}}" if np.isfinite(x) else "N/A"


# ============================================================================
# Core plotting
# ============================================================================


def _make_summary_figure(
    positions: NDArray[np.floating],
    dt: float,
    tau_seconds: NDArray[np.floating],
    msd: NDArray[np.floating],
    alpha: float,
    fit_line_full: NDArray[np.floating],
    d_alpha: NDArray[np.floating],
    tau_vacf_seconds: NDArray[np.floating],
    vacf_norm: NDArray[np.floating],
    tau_p: float,
    tau_p_fit: NDArray[np.floating],
    vacf_fit_curve: NDArray[np.floating],
    tau_orient_seconds: NDArray[np.floating],
    orient_corr: NDArray[np.floating],
    tau_r: float,
    tau_r_fit: NDArray[np.floating],
    orient_fit_curve: NDArray[np.floating],
    freq_psd: NDArray[np.floating],
    psd_vx: NDArray[np.floating],
    psd_vy: NDArray[np.floating],
    beta_x: float,
    beta_y: float,
    fit_freqs_x: NDArray[np.floating],
    fit_line_x: NDArray[np.floating],
    fit_freqs_y: NDArray[np.floating],
    fit_line_y: NDArray[np.floating],
    delta_theta: NDArray[np.floating],
    droplet_id: int | str,
    experiment_name: str,
) -> plt.Figure:
    """
    Assemble the 2 × 3 publication-quality summary figure.

    Parameters
    ----------
    (All arrays are pre-computed by ``analyze_single_droplet``.)
    droplet_id : int or str
        Droplet identifier for titles.
    experiment_name : str
        Dataset label for the overall figure title.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The complete summary figure (not yet saved).
    """
    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    fig.suptitle(
        f"{experiment_name}  —  Droplet {droplet_id}  |  Single-Trajectory Analysis",
        fontsize=14,
        fontweight="bold",
    )
    gs = GridSpec(2, 3, figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])

    # ── (1) Trajectory ───────────────────────────────────────────────────────
    ax1.plot(positions[:, 0], positions[:, 1], lw=1.0, color=StyleTokens.TRAJECTORY, alpha=0.85)
    ax1.plot(*positions[0], "o", ms=7, color=StyleTokens.START, label="Start", zorder=5)
    ax1.plot(*positions[-1], "s", ms=7, color=StyleTokens.END, label="End", zorder=5)
    ax1.set_title("Trajectory")
    ax1.set_xlabel("x (px)")
    ax1.set_ylabel("y (px)")
    ax1.set_aspect("equal")
    ax1.legend(loc="best")
    apply_grid(ax1, alpha=0.25)

    # ── (2) MSD log-log ───────────────────────────────────────────────────────
    ax2.loglog(tau_seconds, msd, color=StyleTokens.MSD, label="MSD")
    fit_mask = ~np.isnan(fit_line_full)
    if np.any(fit_mask):
        ax2.loglog(
            tau_seconds[fit_mask],
            fit_line_full[fit_mask],
            "k:",
            lw=2.0,
            label=rf"$\tau^{{{_safe_label(alpha)}}}$",
        )
    ax2.set_title(rf"MSD  ($\alpha$ = {_safe_label(alpha)})")
    ax2.set_xlabel(r"Lag time $\tau$ (s)")
    ax2.set_ylabel(r"MSD (px²)")
    ax2.legend()
    apply_grid(ax2, alpha=0.25, which="both")

    # ── (3) Local MSD exponent ────────────────────────────────────────────────
    valid = slice(1, -1)
    ax3.semilogx(
        tau_seconds[valid], d_alpha[valid],
        color=StyleTokens.LOCAL_ALPHA, alpha=0.7, label=r"$\alpha(\tau)$"
    )
    ax3.axhline(alpha, ls="--", color="k", lw=1.3,
                label=rf"Global $\alpha$ = {_safe_label(alpha)}")
    ax3.axhline(1.0, ls=":", color="gray", lw=1.0, label="Diffusive (α=1)")
    ax3.axhline(2.0, ls=":", color="gray", lw=1.0, alpha=0.5, label="Ballistic (α=2)")
    ax3.set_title(r"Local MSD Exponent $\alpha(\tau)$")
    ax3.set_xlabel(r"Lag time $\tau$ (s)")
    ax3.set_ylabel(r"$\alpha(\tau)$")
    ax3.legend(fontsize=8)
    apply_grid(ax3, alpha=0.25, which="both")

    # ── (4) VACF ──────────────────────────────────────────────────────────────
    cut = max(len(tau_vacf_seconds) // 4, 1)
    ax4.plot(tau_vacf_seconds[:cut], vacf_norm[:cut],
             color=StyleTokens.VACF, label=r"$C_v(\tau)/C_v(0)$")
    if len(tau_p_fit) > 0 and np.isfinite(tau_p):
        mask = tau_p_fit <= tau_vacf_seconds[cut - 1]
        if np.any(mask):
            ax4.plot(
                tau_p_fit[mask], vacf_fit_curve[mask],
                "k--", lw=1.5,
                label=rf"$\tau_p$ = {_safe_label(tau_p)} s",
            )
    ax4.set_title(rf"VACF  ($\tau_p$ = {_safe_label(tau_p)} s)")
    ax4.set_xlabel(r"Lag time $\tau$ (s)")
    ax4.set_ylabel(r"$C_v(\tau)/C_v(0)$")
    ax4.legend()
    apply_grid(ax4, alpha=0.25)

    # ── (5) Orientation correlation ───────────────────────────────────────────
    ax5.plot(tau_orient_seconds, orient_corr,
             color=StyleTokens.ORIENT_CORR, label=r"$\langle\cos\Delta\theta\rangle$")
    if len(tau_r_fit) > 0 and np.isfinite(tau_r):
        ax5.plot(
            tau_r_fit, orient_fit_curve,
            "k--", lw=1.5,
            label=rf"$\tau_r$ = {_safe_label(tau_r)} s",
        )
    ax5.set_title(rf"Orientation Corr.  ($\tau_r$ = {_safe_label(tau_r)} s)")
    ax5.set_xlabel(r"Lag time $\tau$ (s)")
    ax5.set_ylabel(r"$C_\theta(\tau)$")
    ax5.legend()
    apply_grid(ax5, alpha=0.25)

    # ── (6) PSD vx and vy ────────────────────────────────────────────────────
    if len(freq_psd) > 0:
        ax6.loglog(freq_psd, psd_vx, color=StyleTokens.PSD_X, alpha=0.8, label=r"$S_{v_x}$")
        ax6.loglog(freq_psd, psd_vy, color=StyleTokens.PSD_Y, alpha=0.8, label=r"$S_{v_y}$")
        if len(fit_freqs_x) > 0 and np.isfinite(beta_x):
            ax6.loglog(fit_freqs_x, fit_line_x, "b:", lw=1.8,
                       label=rf"$f^{{-{_safe_label(beta_x, '.2f')}}}$ (x)")
        if len(fit_freqs_y) > 0 and np.isfinite(beta_y):
            ax6.loglog(fit_freqs_y, fit_line_y, "r:", lw=1.8,
                       label=rf"$f^{{-{_safe_label(beta_y, '.2f')}}}$ (y)")
        ax6.set_xlabel("Frequency (Hz)")
        ax6.set_ylabel(r"$S_v(f)$")
    else:
        ax6.text(0.5, 0.5, "PSD: insufficient data",
                 ha="center", va="center", transform=ax6.transAxes, color="gray")
    ax6.set_title(
        rf"PSD  ($\beta_x$ = {_safe_label(beta_x, '.2f')},"
        rf" $\beta_y$ = {_safe_label(beta_y, '.2f')})"
    )
    ax6.legend(fontsize=8)
    apply_grid(ax6, alpha=0.25, which="both")

    return fig


# ============================================================================
# Public API
# ============================================================================


def analyze_single_droplet(
    positions: NDArray[np.float64],
    dt: float,
    droplet_id: int | str,
    experiment_name: str,
    output_dir: Path,
    msd_lag_fraction: float = 0.4,
    msd_fit_fraction: float = 0.3,
    vacf_max_fit_fraction: Optional[float] = None,
    orient_max_fit_fraction: Optional[float] = None,
    psd_fmin: float = 0.1,
    psd_fmax: float = 2.0,
    save_figure: bool = True,
    dpi: int = 300,
) -> SingleDropletAnalysisResult:
    """
    Characterise a single droplet trajectory with publication-quality output.

    All computations delegate to the canonical ``analysis/`` package
    (``analysis.msd``, ``analysis.vacf``, ``analysis.orientation``,
    ``analysis.psd``).  No quantities are recomputed or reimplemented here.

    Parameters
    ----------
    positions : ndarray, shape (N, 2)
        Ordered ``(x, y)`` positions.  Must have N ≥ 8.
    dt : float
        Time step between consecutive frames (seconds).
    droplet_id : int or str
        Identifier used in figure titles and the output filename.
    experiment_name : str
        Dataset label used in figure titles and the output filename.
    output_dir : Path
        Directory where the summary figure is saved.  Created if absent.
    msd_lag_fraction : float
        Fraction of N-1 used as the maximum MSD lag.
        Default ``0.4`` retains 40 % of the lag range.
    msd_fit_fraction : float
        Fraction of the MSD curve included in the log-log power-law fit.
        Default ``0.3`` fits the first 30 % of the lag range, avoiding
        the noisy long-lag tail.
    vacf_max_fit_fraction : float or None
        Maximum fraction of the VACF used for the exponential persistence
        fit.  ``None`` lets the fit stop at the first zero-crossing.
    orient_max_fit_fraction : float or None
        Maximum fraction of the orientation correlation used for the
        exponential fit.
    psd_fmin : float
        Lower frequency bound for the PSD power-law fit (Hz).
    psd_fmax : float
        Upper frequency bound for the PSD power-law fit (Hz).
    save_figure : bool
        If True, write the summary figure to ``output_dir``.
    dpi : int
        Resolution of the saved figure.

    Returns
    -------
    result : SingleDropletAnalysisResult
        Dataclass containing all scalar results and the figure path.

    Raises
    ------
    ValueError
        If ``positions.ndim != 2`` or ``positions.shape[1] != 2``.
        If ``dt <= 0``.
        If ``N < 8`` (too few frames for meaningful analysis).

    Notes
    -----
    **Velocity** is computed once via ``analysis.vacf.compute_velocity``
    and reused for VACF, heading angles, and PSD.  This is the same
    forward-difference definition used everywhere in the ensemble pipeline.

    **Time convention:** All lag arrays returned by the physics functions
    are in frame units.  This function converts them to seconds before
    fitting and plotting (``tau_seconds = tau * dt``), matching the
    convention in ``analysis/main.py``.

    **MSD consistency:** ``fit_msd_power_law`` is called on the same MSD
    array inside ``analyze_single_droplet`` — the identical function used
    by the ensemble pipeline — so the single-droplet ``alpha`` is
    guaranteed to be consistent with ensemble batch results.
    """
    # ── Validate inputs ───────────────────────────────────────────────────────
    positions = np.asarray(positions, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError(
            f"positions must have shape (N, 2); got {positions.shape}"
        )
    n_frames = positions.shape[0]
    if n_frames < 8:
        raise ValueError(
            f"Trajectory must have at least 8 frames; got {n_frames}."
        )
    if dt <= 0:
        raise ValueError(f"dt must be positive; got {dt}.")

    output_dir = Path(output_dir)

    # ── Velocity ──────────────────────────────────────────────────────────────
    # canonical: v = (r[1:] - r[:-1]) / dt   shape (N-1, 2)
    v = compute_velocity(positions, dt)
    speed = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2)
    mean_speed = float(np.mean(speed))
    rms_speed = float(np.sqrt(np.mean(speed ** 2)))

    # ── Heading angles and turning angles ─────────────────────────────────────
    theta = compute_heading_angle(positions, dt)       # (N-1,) unwrapped radians
    delta_theta = _turning_angles(theta)               # (N-2,)
    turn_mean = float(np.mean(delta_theta))
    turn_std = float(np.std(delta_theta, ddof=1)) if len(delta_theta) > 1 else float(np.nan)

    # ── MSD ───────────────────────────────────────────────────────────────────
    msd = compute_msd(positions, lag_fraction=msd_lag_fraction)
    alpha, alpha_err, r2_msd, fit_line_full = fit_msd_power_law(
        msd, fit_fraction=msd_fit_fraction
    )
    tau_msd = np.arange(1, len(msd) + 1)
    tau_seconds = tau_msd * dt

    # ── Local MSD exponent ────────────────────────────────────────────────────
    d_alpha = compute_local_msd_exponent(tau_seconds, msd)

    # ── VACF ──────────────────────────────────────────────────────────────────
    vacf_raw = compute_vacf(positions, dt)             # (N-1,)
    tau_vacf = np.arange(len(vacf_raw))
    tau_vacf_seconds = tau_vacf * dt

    vacf_norm = (
        vacf_raw / vacf_raw[0]
        if vacf_raw[0] != 0
        else vacf_raw.copy()
    )

    tau_p, tau_p_err, tau_p_fit, vacf_fit_curve = fit_vacf_exponential(
        tau_vacf_seconds, vacf_norm,
        max_fit_fraction=vacf_max_fit_fraction,
    )

    # ── Orientation correlation ───────────────────────────────────────────────
    orient_corr = compute_orientation_corr(theta)      # (N-1,)
    tau_orient = np.arange(len(orient_corr))
    tau_orient_seconds = tau_orient * dt

    tau_r, tau_r_err, tau_r_fit, orient_fit_curve = fit_orientation_persistence(
        tau_orient_seconds, orient_corr,
        max_fit_fraction=orient_max_fit_fraction,
    )

    # ── PSD ───────────────────────────────────────────────────────────────────
    freq_psd, psd_vx, psd_vy = compute_psd_components(positions, dt)

    if len(freq_psd) > 0:
        beta_x, bx_err, r2_bx, fit_freqs_x, fit_line_x = fit_psd_powerlaw(
            freq_psd, psd_vx, fmin=psd_fmin, fmax=psd_fmax
        )
        beta_y, by_err, r2_by, fit_freqs_y, fit_line_y = fit_psd_powerlaw(
            freq_psd, psd_vy, fmin=psd_fmin, fmax=psd_fmax
        )
    else:
        beta_x = beta_y = float(np.nan)
        fit_freqs_x = fit_freqs_y = np.array([])
        fit_line_x = fit_line_y = np.array([])

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = _make_summary_figure(
        positions=positions,
        dt=dt,
        tau_seconds=tau_seconds,
        msd=msd,
        alpha=alpha,
        fit_line_full=fit_line_full,
        d_alpha=d_alpha,
        tau_vacf_seconds=tau_vacf_seconds,
        vacf_norm=vacf_norm,
        tau_p=tau_p,
        tau_p_fit=tau_p_fit,
        vacf_fit_curve=vacf_fit_curve,
        tau_orient_seconds=tau_orient_seconds,
        orient_corr=orient_corr,
        tau_r=tau_r,
        tau_r_fit=tau_r_fit,
        orient_fit_curve=orient_fit_curve,
        freq_psd=freq_psd,
        psd_vx=psd_vx,
        psd_vy=psd_vy,
        beta_x=beta_x,
        beta_y=beta_y,
        fit_freqs_x=fit_freqs_x,
        fit_line_x=fit_line_x,
        fit_freqs_y=fit_freqs_y,
        fit_line_y=fit_line_y,
        delta_theta=delta_theta,
        droplet_id=droplet_id,
        experiment_name=experiment_name,
    )

    figure_path: Optional[Path] = None
    if save_figure:
        output_dir.mkdir(parents=True, exist_ok=True)
        fname = f"{experiment_name}_droplet_{droplet_id}_summary.png"
        figure_path = output_dir / fname
        fig.savefig(figure_path, dpi=dpi, bbox_inches="tight")
        print(f"Saved : {figure_path}")

    plt.close(fig)

    # ── Result ────────────────────────────────────────────────────────────────
    return SingleDropletAnalysisResult(
        droplet_id=droplet_id,
        experiment_name=experiment_name,
        track_length=n_frames,
        mean_speed=mean_speed,
        rms_speed=rms_speed,
        msd_alpha=float(alpha),
        msd_alpha_error=float(alpha_err),
        msd_r2=float(r2_msd),
        persistence_time=float(tau_p),
        rotational_persistence_time=float(tau_r),
        psd_beta_x=float(beta_x),
        psd_beta_y=float(beta_y),
        turning_angle_mean=turn_mean,
        turning_angle_std=turn_std,
        figure_path=figure_path,
    )


# ============================================================================
# Analyze All Droplets
# ============================================================================

def analyze_all_droplets(
    positions_dict: Dict[int, NDArray[np.float64]],
    dt: float,
    experiment_name: str,
    output_dir: Path,
    msd_lag_fraction: float = 0.4,
    msd_fit_fraction: float = 0.3,
    vacf_max_fit_fraction: Optional[float] = None,
    orient_max_fit_fraction: Optional[float] = None,
    psd_fmin: float = 0.1,
    psd_fmax: float = 2.0,
    dpi: int = 300,
) -> pd.DataFrame:
    """
    Compare all retained droplets from one experiment.
    
    Produces 7 comparison figures (Trajectory, MSD, Local alpha, VACF, 
    Orientation Correlation, PSD_x, PSD_y) with one curve per droplet.
    Saves a summary statistics CSV.
    
    Parameters
    ----------
    positions_dict : Dict[int, ndarray]
        Mapping from droplet ID to (N, 2) position array.
    dt : float
        Time step between frames in seconds.
    experiment_name : str
        Label for titles and filenames.
    output_dir : Path
        Directory to save plots and statistics CSV.
    msd_lag_fraction : float
        Fraction of N-1 used as the maximum MSD lag.
    msd_fit_fraction : float
        Fraction of the MSD curve included in the log-log power-law fit.
    vacf_max_fit_fraction : float or None
        Maximum fraction of the VACF used for the exponential fit.
    orient_max_fit_fraction : float or None
        Maximum fraction of the orientation correlation used for the exponential fit.
    psd_fmin : float
        Lower frequency bound for PSD fit (Hz).
    psd_fmax : float
        Upper frequency bound for PSD fit (Hz).
    dpi : int
        Figure resolution.
        
    Returns
    -------
    df : pd.DataFrame
        Table of statistics with one row per droplet.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    records = []
    
    # Dictionaries to store curves for plotting
    curves_msd = {}
    curves_dalpha = {}
    curves_vacf = {}
    curves_orient = {}
    curves_psdx = {}
    curves_psdy = {}
    
    for did, pos in positions_dict.items():
        pos = np.asarray(pos, dtype=np.float64)
        n_frames = pos.shape[0]
        if n_frames < 8:
            continue
            
        # Velocity
        v = compute_velocity(pos, dt)
        speed = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2)
        mean_speed = float(np.mean(speed))
        rms_speed = float(np.sqrt(np.mean(speed ** 2)))
        
        # MSD
        msd = compute_msd(pos, lag_fraction=msd_lag_fraction)
        alpha, alpha_err, r2_msd, fit_line_full = fit_msd_power_law(msd, fit_fraction=msd_fit_fraction)
        tau_msd = np.arange(1, len(msd) + 1)
        tau_seconds = tau_msd * dt
        curves_msd[did] = (tau_seconds, msd, alpha)
        
        # Local alpha
        d_alpha = compute_local_msd_exponent(tau_seconds, msd)
        curves_dalpha[did] = (tau_seconds, d_alpha, alpha)
        
        # VACF
        vacf_raw = compute_vacf(pos, dt)
        tau_vacf = np.arange(len(vacf_raw))
        tau_vacf_seconds = tau_vacf * dt
        vacf_norm = vacf_raw / vacf_raw[0] if vacf_raw[0] != 0 else vacf_raw.copy()
        tau_p, tau_p_err, tau_p_fit, vacf_fit_curve = fit_vacf_exponential(tau_vacf_seconds, vacf_norm, max_fit_fraction=vacf_max_fit_fraction)
        curves_vacf[did] = (tau_vacf_seconds, vacf_norm, tau_p)
        
        # Orientation
        theta = compute_heading_angle(pos, dt)
        orient_corr = compute_orientation_corr(theta)
        tau_orient = np.arange(len(orient_corr))
        tau_orient_seconds = tau_orient * dt
        tau_r, tau_r_err, tau_r_fit, orient_fit_curve = fit_orientation_persistence(tau_orient_seconds, orient_corr, max_fit_fraction=orient_max_fit_fraction)
        curves_orient[did] = (tau_orient_seconds, orient_corr, tau_r)
        
        # PSD
        freq_psd, psd_vx, psd_vy = compute_psd_components(pos, dt)
        if len(freq_psd) > 0:
            beta_x, bx_err, r2_bx, fit_freqs_x, fit_line_x = fit_psd_powerlaw(freq_psd, psd_vx, fmin=psd_fmin, fmax=psd_fmax)
            beta_y, by_err, r2_by, fit_freqs_y, fit_line_y = fit_psd_powerlaw(freq_psd, psd_vy, fmin=psd_fmin, fmax=psd_fmax)
        else:
            beta_x = beta_y = float(np.nan)
        curves_psdx[did] = (freq_psd, psd_vx, beta_x)
        curves_psdy[did] = (freq_psd, psd_vy, beta_y)
        
        # Save record
        records.append({
            "droplet_id": did,
            "track_length": n_frames,
            "mean_speed": mean_speed,
            "rms_speed": rms_speed,
            "alpha": float(alpha),
            "alpha_error": float(alpha_err),
            "alpha_r2": float(r2_msd),
            "tau_p": float(tau_p),
            "tau_r": float(tau_r),
            "beta_x": float(beta_x),
            "beta_y": float(beta_y),
        })
        
    df = pd.DataFrame(records)
    csv_path = output_dir / f"{experiment_name}_droplet_statistics.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved : {csv_path}")
    
    # ── Figure 1: Trajectory ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    for did, pos in positions_dict.items():
        if len(pos) < 8: continue
        color = get_droplet_color(did)
        ax.plot(pos[:, 0], pos[:, 1], lw=1.5, color=color, label=f"Droplet {did}")
        ax.plot(*pos[0], "o", ms=6, color=color)
        ax.plot(*pos[-1], "s", ms=6, color=color)
    ax.set_title(f"{experiment_name}\nTrajectory Comparison")
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    ax.set_aspect("equal")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    apply_grid(ax, alpha=0.25)
    fig.savefig(output_dir / "trajectory_comparison.png", dpi=dpi, bbox_inches="tight")
    print(f"Saved : {output_dir / 'trajectory_comparison.png'}")
    plt.close(fig)

    # ── Figure 2: MSD Comparison ───────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    for did, (t_sec, msd, a) in curves_msd.items():
        ax.loglog(t_sec, msd, color=get_droplet_color(did), label=f"Droplet {did} (α={_safe_label(a, '.2f')})")
    ax.set_title("MSD Comparison")
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"MSD (px²)")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    apply_grid(ax, alpha=0.25, which="both")
    fig.savefig(output_dir / "msd_comparison.png", dpi=dpi, bbox_inches="tight")
    print(f"Saved : {output_dir / 'msd_comparison.png'}")
    plt.close(fig)

    # ── Figure 3: Local MSD Exponent Comparison ────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    for did, (t_sec, d_a, a) in curves_dalpha.items():
        valid = slice(1, -1)
        ax.semilogx(t_sec[valid], d_a[valid], color=get_droplet_color(did), label=f"Droplet {did} (α={_safe_label(a, '.2f')})")
    ax.axhline(1.0, ls="--", color="gray", lw=1.5, label="Diffusive (α=1)")
    ax.axhline(2.0, ls=":", color="gray", lw=1.5, label="Ballistic (α=2)")
    ax.set_title("Local MSD Exponent Comparison")
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"$\alpha(\tau)$")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    apply_grid(ax, alpha=0.25, which="both")
    fig.savefig(output_dir / "local_alpha_comparison.png", dpi=dpi, bbox_inches="tight")
    print(f"Saved : {output_dir / 'local_alpha_comparison.png'}")
    plt.close(fig)
    
    # ── Figure 4: VACF Comparison ──────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    for did, (t_sec, vacf_n, tp) in curves_vacf.items():
        cut = max(len(t_sec) // 4, 1)
        ax.plot(t_sec[:cut], vacf_n[:cut], color=get_droplet_color(did), label=f"Droplet {did} (τp={_safe_label(tp, '.1f')} s)")
    ax.axhline(0, ls="--", color="k", lw=1.0)
    ax.set_title("VACF Comparison")
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"$C_v(\tau)/C_v(0)$")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    apply_grid(ax, alpha=0.25)
    fig.savefig(output_dir / "vacf_comparison.png", dpi=dpi, bbox_inches="tight")
    print(f"Saved : {output_dir / 'vacf_comparison.png'}")
    plt.close(fig)

    # ── Figure 5: Orientation Correlation Comparison ───────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    for did, (t_sec, or_corr, tr) in curves_orient.items():
        ax.plot(t_sec, or_corr, color=get_droplet_color(did), label=f"Droplet {did} (τr={_safe_label(tr, '.1f')} s)")
    ax.axhline(0, ls="--", color="k", lw=1.0)
    ax.set_title("Orientation Correlation Comparison")
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"$C_\theta(\tau)$")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    apply_grid(ax, alpha=0.25)
    fig.savefig(output_dir / "orientation_comparison.png", dpi=dpi, bbox_inches="tight")
    print(f"Saved : {output_dir / 'orientation_comparison.png'}")
    plt.close(fig)

    # ── Figure 6: PSD(vx) Comparison ───────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    for did, (f_psd, p_x, bx) in curves_psdx.items():
        if len(f_psd) > 0:
            ax.loglog(f_psd, p_x, color=get_droplet_color(did), label=f"Droplet {did} (βx={_safe_label(bx, '.2f')})")
    ax.set_title("PSD($v_x$) Comparison")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel(r"$S_{v_x}(f)$")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    apply_grid(ax, alpha=0.25, which="both")
    fig.savefig(output_dir / "psd_vx_comparison.png", dpi=dpi, bbox_inches="tight")
    print(f"Saved : {output_dir / 'psd_vx_comparison.png'}")
    plt.close(fig)

    # ── Figure 7: PSD(vy) Comparison ───────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    for did, (f_psd, p_y, by) in curves_psdy.items():
        if len(f_psd) > 0:
            ax.loglog(f_psd, p_y, color=get_droplet_color(did), label=f"Droplet {did} (βy={_safe_label(by, '.2f')})")
    ax.set_title("PSD($v_y$) Comparison")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel(r"$S_{v_y}(f)$")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    apply_grid(ax, alpha=0.25, which="both")
    fig.savefig(output_dir / "psd_vy_comparison.png", dpi=dpi, bbox_inches="tight")
    print(f"Saved : {output_dir / 'psd_vy_comparison.png'}")
    plt.close(fig)
    
    return df


# ============================================================================
# CLI demo
# ============================================================================

def _demo() -> None:
    """
    Run a quick self-contained demonstration when this module is executed
    directly (``python single_droplet_analysis.py``).

    Generates a synthetic active-particle trajectory (persistent random walk)
    and runs ``analyze_single_droplet`` on it.
    """
    rng = np.random.default_rng(0)

    # Synthetic persistent random walk: exponentially correlated speed
    N = 600
    dt = 0.2          # 5 fps equivalent
    tau_p_true = 4.0  # seconds

    vx = np.zeros(N)
    vy = np.zeros(N)
    noise_amp = np.sqrt(2 * dt / tau_p_true)
    for i in range(1, N):
        vx[i] = vx[i - 1] * np.exp(-dt / tau_p_true) + noise_amp * rng.standard_normal()
        vy[i] = vy[i - 1] * np.exp(-dt / tau_p_true) + noise_amp * rng.standard_normal()

    # Integrate to get positions
    x = np.cumsum(vx) * dt
    y = np.cumsum(vy) * dt
    positions = np.column_stack([x, y])

    result = analyze_single_droplet(
        positions=positions,
        dt=dt,
        droplet_id="demo",
        experiment_name="synthetic_PRW",
        output_dir=Path("output/plots/single_droplet"),
        msd_lag_fraction=0.4,
        msd_fit_fraction=0.3,
        psd_fmin=0.01,
        psd_fmax=1.0,
    )
    print(result)

# ============================================================================
#  Main
# ============================================================================



def main() -> None:
    """
    Run single-droplet analysis for one selected trajectory.
    """

    CSV_PATH = Path("./output_droplet_tracking/csv_files")

    file_path = CSV_PATH / "for_2.1_droplets.csv"

    positions_dict, total_paths, consistent_paths = (
        load_positions(file_path)
    )

    # --------------------------------------------------
    # Select droplet ID
    # --------------------------------------------------

    droplet_id = 2

    positions = positions_dict[droplet_id]

    print(f"--- Analyzing Single Droplet ({droplet_id}) ---")
    result = analyze_single_droplet(
        positions=positions,
        dt=0.2,
        droplet_id=droplet_id,
        experiment_name=file_path.stem,
        output_dir=Path(
            "./output_droplet_tracking/plots/single_droplet"
        ),
    )

    print(result)

    # --------------------------------------------------
    # Compare all retained droplets
    # --------------------------------------------------
    print("\n--- Analyzing All Droplets ---")
    df = analyze_all_droplets(
        positions_dict=positions_dict,
        dt=0.2,
        experiment_name=file_path.stem,
        output_dir=Path("./output_droplet_tracking/plots/single_droplet"),
    )
    print("\nDroplet Statistics DataFrame:")
    print(df.to_string())



if __name__ == "__main__":
    main()

