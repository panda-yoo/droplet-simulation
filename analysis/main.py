"""
analysis.main
=============
Post-processing orchestration: load CSV data, compute quantities,
generate all diagnostic plots, and save results.

Contents
--------
- post_processing      : analyse one folder of tracking output
- run_post_processing  : iterate over subfolders and run post_processing
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ── local imports ────────────────────────────────────────
from analysis.extract import load_positions
from analysis.compute import (
    ensemble_msd_analysis,
    ensemble_speed_analysis,
    orientation_analysis,
    ensemble_vacf_analysis,
    ensemble_orientation_corr_analysis,
    ensemble_psd_analysis,
    fit_psd_powerlaw,
)


# =========================================================
# Post-process a single folder
# =========================================================

def post_processing(
    folder_path: str,
    output_dir: str,
    dt: float = 1.0
) -> None:
    """
    Given a folder containing one or more tracking CSVs, compute
    all diagnostic quantities and save the corresponding plots
    for **each** CSV found.

    Plots generated (per CSV)
    -------------------------
    - Consistent trajectories (x-y paths)
    - Log-log MSD
    - Velocity autocorrelation (normalised, first quarter)
    - Orientation correlation
    - PSD of vx  (with power-law fit)
    - PSD of vy  (with power-law fit)
    - <name>_results.txt with exponent values and path counts
    """

    csv_files = [
        f for f in os.listdir(folder_path)
        if f.endswith(".csv")
    ]

    if len(csv_files) == 0:

        print(f"No CSV found in {folder_path}")

        return

    # ── process each CSV in the folder ───────────────────
    for csv_file in csv_files:

        csv_path = os.path.join(
            folder_path,
            csv_file
        )

        # Use the CSV filename (without extension) as label
        label = os.path.splitext(csv_file)[0]

        print(f"Processing : {label}")

        # ── extract positions from CSV ───────────────────
        data, total_paths, consistent_paths = load_positions(csv_path)

        if len(data) == 0:

            print(f"  No valid trajectories in {label}")

            continue

        # Number of droplets used for ensemble averaging
        n_ensemble = consistent_paths

        print(
            f"  Total paths found    : {total_paths}\n"
            f"  Consistent paths     : {consistent_paths}\n"
            f"  Ensemble average (N) : {n_ensemble}"
        )

        # ── compute all quantities ───────────────────────
        tau, ensemble_msd, alpha = (
            ensemble_msd_analysis(data)
        )

        print("MSD exponent alpha =", alpha)

        t, ensemble_speed, ensemble_v = (
            ensemble_speed_analysis(data, dt)
        )

        theta = orientation_analysis(ensemble_v)

        tau_vacf, vacf_raw, vacf_norm = (
            ensemble_vacf_analysis(data, dt)
        )

        tau_orient, orient_corr = (
            ensemble_orientation_corr_analysis(
                data,
                dt
            )
        )

        freq_psd, ensemble_psd_vx, ensemble_psd_vy = (
            ensemble_psd_analysis(data, dt)
        )

        # =================================================
        # PLOT: Consistent trajectories
        # =================================================

        plt.figure(figsize=(7, 7))

        colors = plt.cm.tab10(
            np.linspace(0, 1, n_ensemble)
        )

        for idx, (droplet_id, pos) in enumerate(data.items()):

            plt.plot(
                pos[:, 0],
                pos[:, 1],
                linewidth=1.0,
                color=colors[idx % len(colors)],
                label=f"Droplet {droplet_id}"
            )

            # Mark start and end
            plt.plot(
                pos[0, 0], pos[0, 1],
                'o', color=colors[idx % len(colors)],
                markersize=6
            )
            plt.plot(
                pos[-1, 0], pos[-1, 1],
                's', color=colors[idx % len(colors)],
                markersize=6
            )

        plt.title(
            f"{label} : Consistent Trajectories\n"
            f"N = {n_ensemble} / {total_paths} total paths"
        )

        plt.xlabel('x')
        plt.ylabel('y')
        plt.gca().set_aspect('equal')
        plt.legend(fontsize=7, loc='best')
        plt.tight_layout()

        plt.savefig(
            os.path.join(
                output_dir,
                f"{label}_trajectories.png"
            ),
            dpi=300
        )

        plt.close()

        # =================================================
        # PLOT: Log-log MSD
        # =================================================

        plt.figure(figsize=(6, 4))

        plt.loglog(tau, ensemble_msd)

        plt.title(
            f"{label} : Log-Log MSD\n"
            f"$\\alpha$ = {alpha:.3f}  |  N = {n_ensemble} droplets"
        )

        plt.xlabel(r'Lag time $\tau$')

        plt.ylabel(
            r'$\langle \mathrm{MSD}(\tau)\rangle$'
        )

        plt.tight_layout()

        plt.savefig(
            os.path.join(
                output_dir,
                f"{label}_loglog_msd.png"
            ),
            dpi=300
        )

        plt.close()

        # =================================================
        # PLOT: Velocity autocorrelation
        # =================================================

        plt.figure(figsize=(7, 4))

        cut = len(tau_vacf) // 4

        plt.plot(
            tau_vacf[:cut],
            vacf_norm[:cut]
        )

        plt.title(
            f"{label} : Velocity Autocorrelation\n"
            f"N = {n_ensemble} droplets"
        )

        plt.xlabel(r'Lag time $\tau$')

        plt.ylabel(r'$C_v(\tau)/C_v(0)$')

        plt.tight_layout()

        plt.savefig(
            os.path.join(
                output_dir,
                f"{label}_vacf.png"
            ),
            dpi=300
        )

        plt.close()

        # =================================================
        # PLOT: Orientation correlation
        # =================================================

        plt.figure(figsize=(7, 4))

        plt.plot(
            tau_orient,
            orient_corr
        )

        plt.title(
            f"{label} : Orientation Correlation\n"
            f"N = {n_ensemble} droplets"
        )

        plt.xlabel(r'Lag time $\tau$')

        plt.ylabel(
            r'$\langle \cos(\Delta \theta)\rangle$'
        )

        plt.tight_layout()

        plt.savefig(
            os.path.join(
                output_dir,
                f"{label}_orientation_corr.png"
            ),
            dpi=300
        )

        plt.close()

        # =================================================
        # PLOT: PSD of vx
        # =================================================

        beta_x, fit_line_x = fit_psd_powerlaw(
            freq_psd,
            ensemble_psd_vx
        )

        plt.figure(figsize=(7, 4))

        plt.loglog(
            freq_psd,
            ensemble_psd_vx,
            label='PSD'
        )

        mask_x = (
            (freq_psd >= 0.1)
            &
            (freq_psd <= 2.0)
        )

        plt.loglog(
            freq_psd[mask_x],
            fit_line_x,
            '--',
            label=rf'$f^{{-{beta_x:.2f}}}$'
        )

        plt.title(
            f"{label} : PSD of $v_x$\n"
            f"N = {n_ensemble} droplets"
        )

        plt.xlabel('Frequency')

        plt.ylabel(r'$S_{v_x}(f)$')

        plt.legend()

        plt.tight_layout()

        plt.savefig(
            os.path.join(
                output_dir,
                f"{label}_psd_vx.png"
            ),
            dpi=300
        )

        plt.close()

        # =================================================
        # PLOT: PSD of vy
        # =================================================

        beta_y, fit_line_y = fit_psd_powerlaw(
            freq_psd,
            ensemble_psd_vy
        )

        plt.figure(figsize=(7, 4))

        plt.loglog(
            freq_psd,
            ensemble_psd_vy,
            label='PSD'
        )

        mask_y = (
            (freq_psd >= 0.1)
            &
            (freq_psd <= 2.0)
        )

        plt.loglog(
            freq_psd[mask_y],
            fit_line_y,
            '--',
            label=rf'$f^{{-{beta_y:.2f}}}$'
        )

        plt.title(
            f"{label} : PSD of $v_y$\n"
            f"N = {n_ensemble} droplets"
        )

        plt.xlabel('Frequency')

        plt.ylabel(r'$S_{v_y}(f)$')

        plt.legend()

        plt.tight_layout()

        plt.savefig(
            os.path.join(
                output_dir,
                f"{label}_psd_vy.png"
            ),
            dpi=300
        )

        plt.close()

        # =================================================
        # SAVE: results txt
        # =================================================

        with open(
            os.path.join(
                output_dir,
                f"{label}_results.txt"
            ),
            "w"
        ) as f:

            f.write(
                f"=== {label} ===\n\n"
            )

            f.write(
                f"Total paths detected       = {total_paths}\n"
            )

            f.write(
                f"Consistent paths (used)    = {consistent_paths}\n"
            )

            f.write(
                f"Ensemble average N         = {n_ensemble}\n\n"
            )

            f.write(
                f"MSD exponent alpha         = {alpha}\n"
            )

            f.write(
                f"PSD exponent beta_x        = {beta_x}\n"
            )

            f.write(
                f"PSD exponent beta_y        = {beta_y}\n"
            )

        print(f"{label} PSD exponent beta_x = {beta_x}\n")
        print(f"{label} PSD exponent beta_y = {beta_y}\n")

        print(
            f"Saved results for : {label}"
        )


# =========================================================
# Iterate over folders and run post-processing
# =========================================================

def run_post_processing(
    base_dir: str,
    output_dir: str,
    dt: float = 1.0
) -> None:
    """
    Find CSVs and run ``post_processing``.

    Handles two layouts:
      1. CSVs directly in *base_dir*  → process base_dir itself
      2. CSVs in sub-folders of *base_dir* → process each sub-folder
    """

    if not os.path.exists(base_dir):
        print(f"Base directory '{base_dir}' does not exist. Skipping post-processing.")
        return

    os.makedirs(output_dir, exist_ok=True)

    # ── Check if CSVs exist directly in base_dir ─────────
    csv_in_base = [
        f for f in os.listdir(base_dir)
        if f.endswith(".csv") and os.path.isfile(
            os.path.join(base_dir, f)
        )
    ]

    if csv_in_base:
        # CSVs are directly here — process base_dir itself
        post_processing(
            base_dir,
            output_dir,
            dt
        )

    # ── Also check sub-folders ───────────────────────────
    for folder in os.listdir(base_dir):

        folder_path = os.path.join(
            base_dir,
            folder
        )

        if os.path.isdir(folder_path):

            post_processing(
                folder_path,
                output_dir,
                dt
            )

