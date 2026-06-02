"""
main.py — Unified Entry Point
==============================
Provides the ``DropletPipeline`` class for running simulation, tracking,
and analysis from a single place.

Typical workflow
----------------
1. Run tracking once to generate CSV files from video:

       pipeline.run_tracking()

2. Run analysis repeatedly on those CSV files:

       pipeline.run_analysis(base_dir="output_droplet_tracking/csv_files",
                             output_dir="output_droplet_tracking/plots",
                             dt=0.2)

3. Optionally run the simulation to generate synthetic trajectories,
   then analyse them:

       pipeline.run_simulation()
       pipeline.run_analysis(base_dir="simulation data",
                             output_dir="./output_sim",
                             dt=0.002)

Most users should call ``pipeline.run_analysis()`` on the command line or
in a notebook.  Simulation and tracking are run much less frequently and
are kept separate to avoid accidental re-execution.

Both simulation and tracking produce CSV files with columns:
    droplet_id, x, y, vx, vy, theta

These CSVs are fed into ``run_analysis``, which computes MSD, VACF, PSD,
orientation correlation, and ensemble-bias diagnostics using the canonical
functions in the ``analysis/`` package.
"""


# ====================================================
# ANALYSIS
# ====================================================

class DropletPipeline:
    """
    Unified interface for the active-droplet project.

    There are four modes of operation, in order of typical usage:

    run_analysis(base_dir, output_dir, dt, ...)
        **Most commonly used.**
        Loads existing trajectory CSVs from ``base_dir`` and runs the
        full analysis pipeline (MSD, VACF, PSD, orientation correlation,
        persistence fits, ensemble-bias diagnostics).  Saves plots and a
        ``<label>_results.txt`` summary to ``output_dir``.

        Use this after ``run_tracking`` has already produced CSVs, or
        any time you want to re-run or tune the analysis without
        re-tracking.

        Example
        -------
        >>> pipeline = DropletPipeline()
        >>> pipeline.run_analysis(
        ...     base_dir="output_droplet_tracking/csv_files",
        ...     output_dir="output_droplet_tracking/plots",
        ...     dt=0.2,               # seconds per frame (30 fps video)
        ...     msd_lag_fraction=0.3,
        ...     msd_fit_fraction=0.3,
        ... )

    run_tracking(output_dir, save_annotated)
        Run YOLO-based object detection and tracking on every video in
        ``./videos/``.  Produces per-video trajectory CSVs and calls
        ``run_analysis`` automatically on the generated files.

        Run this once per batch of videos.  After the CSVs exist, use
        ``run_analysis`` directly to avoid re-running tracking.

        Example
        -------
        >>> pipeline.run_tracking(save_annotated=True)

    run_simulation()
        Execute the FEniCSx PDE + SDE active-droplet simulation.
        Generates a trajectory CSV at ``simulation data/droplet_<N>_tracking.csv``
        together with an animated GIF and trajectory plots.

        Example
        -------
        >>> pipeline.run_simulation()
        >>> pipeline.run_analysis(base_dir="simulation data",
        ...                       output_dir="./output_sim",
        ...                       dt=0.002)

    run_all(output_dir, save_annotated)
        Convenience wrapper: runs ``run_simulation`` then ``run_tracking``
        (which includes analysis).  Rarely needed outside of a fresh-start
        end-to-end test.

        Example
        -------
        >>> pipeline.run_all()
    """

    # ────────────────────────────────────────────────────────
    # ANALYSIS
    # ────────────────────────────────────────────────────────

    def run_analysis(
        self,
        base_dir: str,
        output_dir: str = "./output",
        dt: float = 0.2,
        msd_lag_fraction: float = 0.3,
        msd_fit_fraction: float = 0.4,
        msd_minimum_pairs: int | None = None,
        vacf_max_fit_fraction: float | None = None,
        orient_max_fit_fraction: float | None = None,
    ):
        """
        Run ONLY the analysis / plotting on existing trajectory CSVs.

        All physics is computed by the canonical ``analysis/`` modules
        (``analysis.msd``, ``analysis.vacf``, ``analysis.orientation``,
        ``analysis.psd``, ``analysis.diagnostics``).  The MSD exponent
        computed here and the exponent shown in plots are guaranteed to be
        numerically identical because both use the same ``fit_msd_power_law``
        call path.

        Parameters
        ----------
        base_dir : str
            Folder (or parent of sub-folders) containing ``.csv`` files.
            Each sub-folder that contains a CSV is processed independently.
        output_dir : str
            Where to save plots and ``<label>_results.txt``.
        dt : float
            Time step between consecutive frames in seconds.
            Must match the frame rate used during tracking or the
            ``dt_val`` used during simulation.
        msd_lag_fraction : float
            Fraction of each trajectory length used as the maximum MSD lag.
        msd_fit_fraction : float
            Fraction of the MSD curve included in the log-log power-law fit.
        msd_minimum_pairs : int or None
            Minimum number of displacement pairs required at every lag.
        vacf_max_fit_fraction : float or None
            Fraction of the VACF used for the exponential persistence fit.
        orient_max_fit_fraction : float or None
            Fraction of the orientation correlation used for the
            exponential persistence fit.

        Example
        -------
        >>> pipeline = DropletPipeline()
        >>> pipeline.run_analysis(
        ...     base_dir="./output_droplet_tracking",
        ...     output_dir="./output_droplet_tracking/plots",
        ...     dt=0.2,
        ... )
        """
        from analysis.main import run_post_processing

        print("=" * 60)
        print("  RUNNING ANALYSIS")
        print("=" * 60)

        run_post_processing(
            base_dir=base_dir,
            output_dir=output_dir,
            dt=dt,
            msd_lag_fraction=msd_lag_fraction,
            msd_fit_fraction=msd_fit_fraction,
            msd_minimum_pairs=msd_minimum_pairs,
            vacf_max_fit_fraction=vacf_max_fit_fraction,
            orient_max_fit_fraction=orient_max_fit_fraction,
        )

    # ────────────────────────────────────────────────────────
    # TRACKING
    # ────────────────────────────────────────────────────────

    def run_tracking(
        self,
        output_dir: str = "output_droplet_tracking",
        save_annotated: bool = True,
    ):
        """
        Execute YOLO-based video tracking and post-processing.

        Discovers all videos in ``./videos/``, runs YOLO object tracking,
        writes per-video trajectory CSVs, and calls ``run_analysis``
        automatically on the generated files.

        Parameters
        ----------
        output_dir : str
            Root directory to save all outputs:
            - ``<output_dir>/csv_files/``         — trajectory CSVs
            - ``<output_dir>/annotated_videos/``  — annotated videos
            - ``<output_dir>/plots/``             — analysis plots

        save_annotated : bool
            If True, saves the YOLO-annotated videos with detection boxes.

        Requires
        --------
        - ``ultralytics`` package
        - ``weight/best.pt``  YOLO model file
        - Video files in ``./videos/`` (supported: .mp4, .mov, .avi, .mkv)
        """
        from tracking.main import run_tracking

        print("=" * 60)
        print("  RUNNING TRACKING + ANALYSIS")
        print("=" * 60)

        run_tracking(
            output_dir=output_dir,
            save_annotated=save_annotated,
        )

    # ────────────────────────────────────────────────────────
    # SIMULATION
    # ────────────────────────────────────────────────────────

    def run_simulation(self):
        """
        Execute the full active-droplet simulation.

        Runs the FEniCSx PDE (chemical trail diffusion) + Euler-Maruyama
        SDE (stochastic droplet motion) for ``Ndrop = 6`` droplets on a
        400 × 400 mesh.

        Outputs (written to ``simulation data/``)
        ------------------------------------------
        - ``multi_droplet_active.gif``              — animated field
        - ``final_field.png``                       — chemical field at T
        - ``droplet_<N>_trajectories.png``          — x-y trajectory plot
        - ``droplet_<N>_tracking.csv``              — trajectory CSV
        - ``data_droplet_<N>_trajectories.npz``     — raw trajectory arrays

        To analyse the simulation output afterwards::

            pipeline.run_analysis(
                base_dir="simulation data",
                output_dir="./output_sim",
                dt=0.002,   # must match dt_val in simulation/main.py
            )

        Requires
        --------
        mpi4py, petsc4py, dolfinx, pyvista, tqdm
        """
        from simulation.main import run_simulation

        print("=" * 60)
        print("  RUNNING SIMULATION")
        print("=" * 60)

        run_simulation()

    # ────────────────────────────────────────────────────────
    # FULL PIPELINE
    # ────────────────────────────────────────────────────────

    def run_all(
        self,
        output_dir: str = "output_droplet_tracking",
        save_annotated: bool = True,
    ):
        """
        Run simulation first, then tracking + analysis.

        Equivalent to calling ``run_simulation()`` followed by
        ``run_tracking(output_dir, save_annotated)``.

        Use this for an end-to-end test or a fresh-start run.
        In most research workflows, simulation and tracking are run
        once and analysis is iterated independently.

        Parameters
        ----------
        output_dir : str
            Root directory for tracking and analysis outputs.
        save_annotated : bool
            If True, saves YOLO-annotated videos.
        """
        self.run_simulation()
        self.run_tracking(output_dir=output_dir, save_annotated=save_annotated)


# ====================================================
# CLI ENTRY POINT
# ====================================================

def main():
    """
    Command-line entry point.

    Edit the sections below to select the desired workflow.
    Only one section should be active (uncommented) at a time.
    """
    pipeline = DropletPipeline()

    # ====================================================
    # ANALYSIS ONLY  ←  DEFAULT
    # ====================================================
    # Loads existing trajectory CSVs and runs the full analysis.
    # Use this after tracking or simulation has already generated CSVs.

    # pipeline.run_analysis(
    #     base_dir="./output_droplet_tracking",
    #     output_dir="./output_droplet_tracking/plots",
    #     dt=0.2,
    # )

    # ====================================================
    # TRACKING
    # ====================================================
    # Run YOLO tracking on all videos in ./videos/.
    # Generates CSVs and runs analysis automatically.
    # Uncomment to use.

    # pipeline.run_tracking(
    #     save_annotated=True,
    # )

    # ====================================================
    # SIMULATION
    # ====================================================
    # Run the FEniCSx PDE + SDE simulation.
    # Then call run_analysis with dt=0.002 to analyse the output.
    # Uncomment to use.

    # pipeline.run_simulation()

    # pipeline.run_analysis(
    #     base_dir="simulation data",
    #     output_dir="./output_sim",
    #     dt=0.2,
    # )
    
    pipeline.run_analysis(
        base_dir="./simulation data/forBeta_50_Sigma_0.01",
        output_dir="./output_sim",
        dt=0.002,
    )
    


    # ====================================================
    # FULL PIPELINE
    # ====================================================
    # Runs simulation then tracking + analysis end-to-end.
    # Uncomment to use.

    # pipeline.run_all(
    #     output_dir="output_droplet_tracking",
    #     save_annotated=True,
    # )


if __name__ == "__main__":
    main()
