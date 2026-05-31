"""
simulation.utils
================
Helper functions for the active-droplet simulation.

Contents
--------
- initial_condition   : zero-field initial condition for the PDE
- update_paper_parameters : time-dependent R, alpha, gamma as per the paper
- save_trajectory_csv : export droplet trajectories to CSV
"""

import numpy as np
import csv
from pathlib import Path
from typing import List


# =========================================================
# Initial condition for the chemical field (zero everywhere)
# =========================================================

def initial_condition(x):
    """Return a zero scalar field matching the mesh coordinates."""
    return np.zeros_like(x[0])


# =========================================================
# Time-dependent physical parameters (paper model)
# =========================================================

def update_paper_parameters(t, R0, Rshrink, Rmin, delta, eta):
    """
    Compute the time-dependent droplet radius, secretion rate,
    and friction coefficient.

    Parameters
    ----------
    t       : float  — current simulation time
    R0      : float  — initial droplet radius
    Rshrink : float  — radius shrinkage rate
    Rmin    : float  — minimum allowed radius
    delta   : float  — chemical interaction length scale
    eta     : float  — fluid viscosity

    Returns
    -------
    Rval     : float  — current radius
    alphaval : float  — secretion rate
    gammaval : float  — Stokes friction coefficient
    """
    Rval = max(R0 - Rshrink * t, Rmin)
    alphaval = 3.0 * (Rval ** 2) / (delta ** 3) * Rshrink
    gammaval = 6.0 * np.pi * eta * Rval
    return Rval, alphaval, gammaval


# =========================================================
# Save trajectories to CSV
# =========================================================

def save_trajectory_csv(
    trajectory_x: List[List[float]],
    trajectory_y: List[List[float]],
    dt_val: float,
    Ndrop: int,
    output_path: str | Path,
) -> None:
    """
    Write per-droplet positions, velocities, and heading angles to CSV.

    Parameters
    ----------
    trajectory_x : list of list
        Per-droplet x-positions over time.
    trajectory_y : list of list
        Per-droplet y-positions over time.
    dt_val : float
        Time step between frames.
    Ndrop : int
        Number of droplets.
    output_path : str or Path
        Destination CSV path.

    Returns
    -------
    None
    """
    rows = []

    for i in range(Ndrop):

        x_traj = trajectory_x[i]
        y_traj = trajectory_y[i]

        for j in range(len(x_traj)):

            x = x_traj[j]
            y = y_traj[j]

            # Velocity and angle (forward difference)
            if j < len(x_traj) - 1:
                vx = (x_traj[j + 1] - x_traj[j]) / dt_val
                vy = (y_traj[j + 1] - y_traj[j]) / dt_val
            elif len(x_traj) > 1:
                vx = (x_traj[j] - x_traj[j - 1]) / dt_val
                vy = (y_traj[j] - y_traj[j - 1]) / dt_val
            else:
                vx = 0.0
                vy = 0.0

            theta = np.arctan2(vy, vx) if (vx != 0.0 or vy != 0.0) else 0.0

            rows.append([
                i,
                x,
                y,
                vx,
                vy,
                theta,
            ])

    with open(output_path, "w", newline="") as f:

        writer = csv.writer(f)

        writer.writerow([
            'droplet id',
            'x',
            'y',
            'vx',
            'vy',
            'theta',
        ])

        writer.writerows(rows)

    print(f"CSV saved to: {output_path}")
