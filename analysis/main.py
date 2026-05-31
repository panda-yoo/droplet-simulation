"""
analysis.main
=============
Post-processing orchestration: load CSV data, compute quantities,
generate all diagnostic plots, and save results.

All physical quantities are computed by the canonical modules:
    analysis.msd         — MSD, fit, local exponent, pair-weighted MSD
    analysis.vacf        — VACF, velocity, persistence fit
    analysis.orientation — orientation correlation, persistence fit
    analysis.psd         — PSD, power-law fit
    analysis.diagnostics — bias diagnostics

Time convention
---------------
Internal lag indices are integers (frame units).
All reporting and plotting uses ``tau_seconds = tau * dt``.

Contents
--------
- post_processing      : analyse one folder of tracking output
- run_post_processing  : iterate over subfolders and run post_processing
"""

import os
import numpy as np
from typing import Optional

from analysis.trajectories import load_positions_with_lengths
from analysis.msd import (
    ensemble_msd_analysis,
    compute_local_msd_exponent,
    compute_pair_weighted_msd,
)
from analysis.vacf import (
    ensemble_vacf_analysis,
    fit_vacf_exponential,
)
from analysis.orientation import (
    ensemble_orientation_corr_analysis,
    fit_orientation_persistence,
)
from analysis.psd import (
    ensemble_psd_analysis,
    fit_psd_powerlaw,
)
from analysis.diagnostics import (
    retained_vs_discarded_lengths,
    retained_vs_discarded_speed_stats,
)
from analysis.plotting import (
    plot_trajectories,
    plot_msd_loglog,
    plot_vacf,
    plot_orientation_corr,
    plot_psd_component,
    plot_local_msd_exponent,
    plot_msd_comparison,
    plot_track_length_hist,
)


# =========================================================
# Post-process a single folder
# =========================================================

def post_processing(
    folder_path: str,
    output_dir: str,
    dt: float = 1.0,
    msd_lag_fraction: float = 0.3,
    msd_fit_fraction: float = 0.3,
    msd_minimum_pairs: Optional[int] = None,
    vacf_max_fit_fraction: Optional[float] = None,
    orient_max_fit_fraction: Optional[float] = None,
) -> None:
    """
    Compute all diagnostic quantities for every CSV in ``folder_path``
    and save plots and a results summary.

    Parameters
    ----------
    folder_path : str
        Directory containing one or more tracking CSVs.
    output_dir : str
        Output directory for plots and ``<label>_results.txt``.
    dt : float
        Time step between consecutive frames (seconds).
    msd_lag_fraction : float
        Fraction of each trajectory length used as the maximum MSD lag.
    msd_fit_fraction : float
        Fraction of the MSD curve included in the log-log power-law fit.
        The same value is used by both ``ensemble_msd_analysis`` (batch)
        and the plotting code so exponents are guaranteed identical.
    msd_minimum_pairs : int or None
        Minimum number of displacement pairs required at every MSD lag.
    vacf_max_fit_fraction : float or None
        Fraction of the VACF used for the exponential persistence fit.
    orient_max_fit_fraction : float or None
        Fraction of the orientation correlation used for the exponential fit.

    Notes
    -----
    Plots generated per CSV
    -----------------------
    - ``<label>_trajectories.png``        — x-y paths
    - ``<label>_loglog_msd.png``          — log-log MSD with fit
    - ``<label>_alpha_vs_tau.png``        — local MSD exponent
    - ``<label>_msd_pair_weighted.png``   — ensemble vs pair-weighted MSD
    - ``<label>_vacf.png``                — normalised VACF with fit
    - ``<label>_orientation_corr.png``    — orientation correlation with fit
    - ``<label>_psd_vx.png``              — PSD of vx with power-law fit
    - ``<label>_psd_vy.png``              — PSD of vy with power-law fit
    - ``<label>_track_length_hist.png``   — retained vs discarded lengths
    - ``<label>_results.txt``             — exponent summary
    """
    csv_files = [f for f in os.listdir(folder_path) if f.endswith(".csv")]

    if not csv_files:
        print(f"No CSV found in {folder_path}")
        return

    for csv_file in csv_files:
        csv_path = os.path.join(folder_path, csv_file)
        label = os.path.splitext(csv_file)[0]
        print(f"Processing : {label}")

        # ── load positions ────────────────────────────────
        (
            data,
            total_paths,
            consistent_paths,
            positions_all,
            lengths,
        ) = load_positions_with_lengths(csv_path)

        if not data:
            print(f"  No valid trajectories in {label}")
            continue

        n_ensemble = consistent_paths
        print(
            f"  Total paths found    : {total_paths}\n"
            f"  Consistent paths     : {consistent_paths}\n"
            f"  Ensemble average (N) : {n_ensemble}"
        )

        os.makedirs(output_dir, exist_ok=True)

        # ── MSD ───────────────────────────────────────────
        # ensemble_msd_analysis calls fit_msd_power_law internally.
        # The same fit_fraction and minimum_pairs are used here and in
        # any notebook that calls ensemble_msd_analysis directly, so
        # batch and plot exponents are numerically identical.
        tau, ensemble_msd, alpha, alpha_err, r2_msd, fit_line_msd = (
            ensemble_msd_analysis(
                data,
                lag_fraction=msd_lag_fraction,
                fit_fraction=msd_fit_fraction,
                minimum_pairs=msd_minimum_pairs,
            )
        )
        tau_seconds = tau * dt
        print(f"  MSD exponent alpha   : {alpha:.4f} ± {alpha_err:.4f}  R²={r2_msd:.4f}")

        # ── local MSD exponent ────────────────────────────
        d_msd = compute_local_msd_exponent(tau_seconds, ensemble_msd)

        # ── pair-weighted MSD (diagnostic) ───────────────
        tau_weighted, msd_weighted = compute_pair_weighted_msd(
            positions_all,
            lag_fraction=msd_lag_fraction,
            minimum_pairs=msd_minimum_pairs,
        )
        tau_weighted_seconds = tau_weighted * dt

        # ── VACF ──────────────────────────────────────────
        tau_vacf, vacf_raw, vacf_norm = ensemble_vacf_analysis(data, dt)
        tau_vacf_seconds = tau_vacf * dt

        tau_p, tau_p_err, r2_vacf, tau_p_fit, vacf_fit_curve = fit_vacf_exponential(
            tau_vacf_seconds,
            vacf_norm,
            max_fit_fraction=vacf_max_fit_fraction,
        )

        # ── orientation correlation ───────────────────────
        tau_orient, orient_corr = ensemble_orientation_corr_analysis(data, dt)
        tau_orient_seconds = tau_orient * dt

        tau_r, tau_r_err, r2_orient, tau_r_fit, orient_fit_curve = fit_orientation_persistence(
            tau_orient_seconds,
            orient_corr,
            max_fit_fraction=orient_max_fit_fraction,
        )

        # ── PSD ───────────────────────────────────────────
        freq_psd, ensemble_psd_vx, ensemble_psd_vy = ensemble_psd_analysis(data, dt)

        beta_x, beta_x_err, r2_x, fit_freqs_x, fit_line_x = fit_psd_powerlaw(
            freq_psd, ensemble_psd_vx
        )
        beta_y, beta_y_err, r2_y, fit_freqs_y, fit_line_y = fit_psd_powerlaw(
            freq_psd, ensemble_psd_vy
        )
        print(f"  PSD exponent beta_x  : {beta_x:.4f}")
        print(f"  PSD exponent beta_y  : {beta_y:.4f}")

        # ── diagnostics ───────────────────────────────────
        retained_lengths, discarded_lengths, n_retained, n_discarded = (
            retained_vs_discarded_lengths(lengths)
        )
        speed_stats = retained_vs_discarded_speed_stats(positions_all, lengths, dt)
        print(
            f"  Retained tracks      : {n_retained}\n"
            f"  Discarded tracks     : {n_discarded}"
        )

        # ── plots ─────────────────────────────────────────
        plot_trajectories(data, label, output_dir, total_paths)

        plot_msd_loglog(
            tau_seconds, ensemble_msd, alpha, alpha_err, r2_msd, fit_line_msd,
            label, output_dir, n_ensemble,
        )

        plot_local_msd_exponent(
            tau_seconds, d_msd, alpha, alpha_err, label, output_dir, n_ensemble,
        )

        if len(tau_weighted) > 0:
            diagnostics_dir = os.path.join(output_dir, "diagnostics")
            os.makedirs(diagnostics_dir, exist_ok=True)
            plot_msd_comparison(
                tau_seconds, ensemble_msd,
                tau_weighted_seconds, msd_weighted,
                label, diagnostics_dir, n_ensemble,
            )

        plot_vacf(
            tau_vacf_seconds, vacf_norm, label, output_dir, n_ensemble,
            tau_p_fit=tau_p_fit, vacf_fit_curve=vacf_fit_curve, tau_p=tau_p,
            tau_p_err=tau_p_err, r2=r2_vacf
        )

        plot_orientation_corr(
            tau_orient_seconds, orient_corr, label, output_dir, n_ensemble,
            tau_r_fit=tau_r_fit, orient_fit_curve=orient_fit_curve, tau_r=tau_r,
            tau_r_err=tau_r_err, r2=r2_orient
        )

        if len(freq_psd) > 0:
            plot_psd_component(
                freq_psd, ensemble_psd_vx,
                beta_x, beta_x_err, r2_x, fit_freqs_x, fit_line_x,
                "vx", label, output_dir, n_ensemble,
            )
            plot_psd_component(
                freq_psd, ensemble_psd_vy,
                beta_y, beta_y_err, r2_y, fit_freqs_y, fit_line_y,
                "vy", label, output_dir, n_ensemble,
            )

        if lengths:
            diagnostics_dir = os.path.join(output_dir, "diagnostics")
            os.makedirs(diagnostics_dir, exist_ok=True)
            plot_track_length_hist(
                retained_lengths, discarded_lengths, label, diagnostics_dir, n_ensemble
            )

        # ── results txt ───────────────────────────────────
        results_path = os.path.join(output_dir, f"{label}_results.txt")
        with open(results_path, "w") as f:
            f.write(f"=== {label} ===\n\n")
            f.write(f"Total paths detected        = {total_paths}\n")
            f.write(f"Retained tracks             = {n_retained}\n")
            f.write(f"Discarded tracks            = {n_discarded}\n")
            f.write(f"Consistent paths (used)     = {consistent_paths}\n")
            f.write(f"Ensemble average N          = {n_ensemble}\n\n")
            f.write(f"MSD exponent alpha          = {alpha} ± {alpha_err}  R2 = {r2_msd}\n")
            f.write(f"PSD exponent beta_x         = {beta_x} ± {beta_x_err}  R2 = {r2_x}\n")
            f.write(f"PSD exponent beta_y         = {beta_y} ± {beta_y_err}  R2 = {r2_y}\n")
            f.write(f"Velocity persistence tau_p  = {tau_p} ± {tau_p_err} s\n")
            f.write(f"Rotational persistence tau_r = {tau_r} ± {tau_r_err} s\n\n")
            f.write("--- Bias diagnostics ---\n")
            f.write(f"Retained mean speed   = {speed_stats['retained_mean_speed']:.6g}\n")
            f.write(f"Retained std speed    = {speed_stats['retained_std_speed']:.6g}\n")
            f.write(f"Discarded mean speed  = {speed_stats['discarded_mean_speed']:.6g}\n")
            f.write(f"Discarded std speed   = {speed_stats['discarded_std_speed']:.6g}\n")

        print(f"  Saved results for : {label}")


# =========================================================
# Iterate over folders and run post-processing
# =========================================================

def run_post_processing(
    base_dir: str,
    output_dir: str,
    dt: float = 1.0,
    msd_lag_fraction: float = 0.3,
    msd_fit_fraction: float = 0.3,
    msd_minimum_pairs: Optional[int] = None,
    vacf_max_fit_fraction: Optional[float] = None,
    orient_max_fit_fraction: Optional[float] = None,
) -> None:
    """
    Find CSVs and run ``post_processing``.

    Handles two layouts:
      1. CSVs directly in ``base_dir``  → process ``base_dir`` itself.
      2. CSVs in sub-folders of ``base_dir`` → process each sub-folder.

    Parameters
    ----------
    base_dir : str
        Base directory containing CSVs or sub-folders with CSVs.
    output_dir : str
        Output directory for plots and results.
    dt : float
        Time step between consecutive frames (seconds).
    msd_lag_fraction : float
        Fraction of each trajectory length used as the maximum MSD lag.
    msd_fit_fraction : float
        Fraction of the MSD curve included in the log-log power-law fit.
    msd_minimum_pairs : int or None
        Minimum number of displacement pairs required at every MSD lag.
    vacf_max_fit_fraction : float or None
        Fraction of the VACF used for the exponential persistence fit.
    orient_max_fit_fraction : float or None
        Fraction of the orientation correlation used for the exponential fit.
    """
    if not os.path.exists(base_dir):
        print(f"Base directory '{base_dir}' does not exist. Skipping.")
        return

    os.makedirs(output_dir, exist_ok=True)

    kwargs = dict(
        dt=dt,
        msd_lag_fraction=msd_lag_fraction,
        msd_fit_fraction=msd_fit_fraction,
        msd_minimum_pairs=msd_minimum_pairs,
        vacf_max_fit_fraction=vacf_max_fit_fraction,
        orient_max_fit_fraction=orient_max_fit_fraction,
    )

    # CSVs directly in base_dir
    csv_in_base = [
        f for f in os.listdir(base_dir)
        if f.endswith(".csv") and os.path.isfile(os.path.join(base_dir, f))
    ]
    if csv_in_base:
        post_processing(base_dir, output_dir, **kwargs)

    # Sub-folders
    for folder in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder)
        if os.path.isdir(folder_path):
            post_processing(folder_path, output_dir, **kwargs)
