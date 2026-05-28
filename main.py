"""
main.py — Unified Entry Point
==============================
Provides the ``DropletPipeline`` class for running
simulation, tracking, and analysis from a single place.

Usage
-----
    from main import DropletPipeline

    pipeline = DropletPipeline()

    # Run only the simulation (FEniCSx active-droplet PDE + SDE)
    pipeline.run_simulation()

    # Run only the YOLO-based video tracking + CSV generation + analysis
    pipeline.run_tracking()

    # Or run everything
    pipeline.run_all()

Both simulation and tracking produce CSV files with columns:
    droplet id, x, y, vx, vy, theta

These CSVs can then be fed into the analysis pipeline to compute
MSD, VACF, PSD, orientation correlation, etc.
"""


# =========================================================
# DropletPipeline — main orchestration class
# =========================================================

class DropletPipeline:
    """
    Unified interface for the active-droplet project.

    Methods
    -------
    run_simulation()
        Run the FEniCSx PDE + SDE simulation.
        Outputs: GIF animation, trajectory plot, .npz, .csv

    run_tracking()
        Run YOLO tracking on videos in ./videos/.
        Outputs: annotated videos, per-video CSVs, analysis plots

    run_analysis(base_dir, output_dir, dt)
        Run ONLY the analysis on existing CSVs (no simulation, no tracking).

    run_all()
        Run simulation first, then tracking + analysis.
    """

    # ─────────────────────────────────────────────────────
    # Simulation
    # ─────────────────────────────────────────────────────

    def run_simulation(self):
        """
        Execute the full active-droplet simulation.

        Requires: mpi4py, petsc4py, dolfinx, pyvista
        Produces:
          - multi_droplet_active.gif
          - final_field.png
          - droplet_<N>_trajectories.png
          - simulation data/droplet_<N>_tracking.csv
          - simulation data/data_droplet_<N>_trajectories.npz
        """
        from simulation.main import run_simulation
        print("=" * 60)
        print("  RUNNING SIMULATION")
        print("=" * 60)
        run_simulation()

    # ─────────────────────────────────────────────────────
    # Tracking + Analysis
    # ─────────────────────────────────────────────────────

    def run_tracking(
        self,
        output_dir: str = "output_droplet_tracking",
        save_annotated: bool = True
    ):
        """
        Execute YOLO-based video tracking and post-processing.

        Parameters
        ----------
        output_dir     : str  — main directory to save all outputs (CSVs, videos, plots)
        save_annotated : bool — whether to save the annotated videos with tracking boxes

        Requires: ultralytics, best.pt model file, videos in ./videos/
        Produces:
          - <output_dir>/annotated_videos/ (if save_annotated is True)
          - <output_dir>/csv_files/
          - <output_dir>/plots/ (MSD, VACF, PSD, orientation plots + results.txt)
        """
        from tracking.main import run_tracking
        print("=" * 60)
        print("  RUNNING TRACKING + ANALYSIS")
        print("=" * 60)
        run_tracking(
            output_dir=output_dir,
            save_annotated=save_annotated
        )

    # ─────────────────────────────────────────────────────
    # Analysis only (on existing CSVs)
    # ─────────────────────────────────────────────────────

    def run_analysis(
        self,
        base_dir: str,
        output_dir: str = "./output",
        dt: float = 1.0 / 30.0,
    ):
        """
        Run ONLY the analysis / plotting on existing CSVs.

        Parameters
        ----------
        base_dir   : str   — folder (or parent of sub-folders) containing CSVs.
                             Each sub-folder with a .csv file will be processed.
        output_dir : str   — where to save plots and results.txt
        dt         : float — time step between frames (default 1/30 fps)

        Example
        -------
        >>> pipeline = DropletPipeline()
        >>> # For simulation output:
        >>> pipeline.run_analysis("simulation data", "./output_sim")
        >>> # For tracking output:
        >>> pipeline.run_analysis("output_droplet_tracking/csv_files", "./output")
        """
        from analysis.main import run_post_processing
        print("=" * 60)
        print("  RUNNING ANALYSIS ONLY")
        print("=" * 60)
        run_post_processing(
            base_dir=base_dir,
            output_dir=output_dir,
            dt=dt,
        )

    # ─────────────────────────────────────────────────────
    # Run everything
    # ─────────────────────────────────────────────────────

    def run_all(self, output_dir: str = "output_droplet_tracking", save_annotated: bool = True):
        """Run simulation, then tracking + analysis."""
        self.run_simulation()
        self.run_tracking(output_dir=output_dir, save_annotated=save_annotated)


# =========================================================
# CLI entry point
# =========================================================

if __name__ == "__main__":
    pipeline = DropletPipeline()
    
    # ── Run the Tracking + Tracking Analysis) ──
    pipeline.run_tracking(save_annotated = True)
    pipeline.run_analysis(base_dir='./output_droplet_tracking',dt=0.2)
    # ── Run the entire pipeline (Simulation + Tracking + Tracking Analysis) ──
    # pipeline.run_all(output_dir="output_droplet_tracking", save_annotated=True)

    # ── Finally, run Analysis on the Simulation Data ─────────────────────────
    # pipeline.run_analysis(
    #     base_dir="simulation data",
    #     output_dir="./output_sim",
    #     dt=1.0,
    # )

