# =========================================================
# simulation package
# ---------------------------------------------------------
# Provides the FEniCSx-based active-droplet simulation.
#   simulation.main.run_simulation()  — full simulation loop
#   simulation.utils                  — helper functions
# =========================================================

from simulation.utils import (
    initial_condition,
    update_paper_parameters,
    save_trajectory_csv,
)
