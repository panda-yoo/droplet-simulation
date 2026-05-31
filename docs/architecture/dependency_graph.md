# Dependency Graph

Full call chains from every entry point to every leaf function.

---

## Entry point: `main.py`

```
main.py  (DropletPipeline)
├── run_simulation()
│    └── simulation.main.run_simulation()
│         ├── simulation.utils.initial_condition()
│         ├── simulation.utils.update_paper_parameters()
│         └── simulation.utils.save_trajectory_csv()
│
├── run_tracking(output_dir, save_annotated)
│    └── tracking.main.run_tracking()
│         ├── tracking.utils.generate_csvs()
│         │    ├── analysis.vacf.compute_velocity()        [via _velocity_for_csv]
│         │    └── tracking.utils.positions_to_rows()
│         └── analysis.main.run_post_processing()          [see below]
│
├── run_analysis(base_dir, output_dir, dt, ...)
│    └── analysis.main.run_post_processing()               [see below]
│
└── run_all()
     ├── run_simulation()
     └── run_tracking()
```

---

## `analysis.main.run_post_processing`

```
analysis.main.run_post_processing(base_dir, output_dir, dt, ...)
└── analysis.main.post_processing(folder_path, output_dir, dt, ...)
     │
     ├── analysis.trajectories.load_positions_with_lengths(csv_path)
     │    └── analysis.trajectories.load_positions(csv_path)
     │
     ├── analysis.msd.ensemble_msd_analysis(data, lag_fraction, fit_fraction, ...)
     │    ├── analysis.msd.compute_msd(r, lag_fraction, minimum_pairs)
     │    └── analysis.msd.fit_msd_power_law(ensemble_msd, fit_fraction, ...)
     │         └── analysis.msd._loglog_ols(log_tau, log_msd)
     │
     ├── analysis.msd.compute_local_msd_exponent(tau_seconds, ensemble_msd)
     │
     ├── analysis.msd.compute_pair_weighted_msd(positions_all, lag_fraction, ...)
     │
     ├── analysis.vacf.ensemble_vacf_analysis(data, dt)
     │    └── analysis.vacf.compute_vacf(r, dt)
     │         └── analysis.vacf.compute_velocity(r, dt)
     │
     ├── analysis.vacf.fit_vacf_exponential(tau_vacf_seconds, vacf_norm, ...)
     │
     ├── analysis.orientation.ensemble_orientation_corr_analysis(data, dt)
     │    ├── analysis.vacf.compute_velocity(r, dt)        [imported into orientation]
     │    ├── analysis.orientation.compute_heading_angle(r, dt)
     │    └── analysis.orientation.compute_orientation_corr(theta)
     │
     ├── analysis.orientation.fit_orientation_persistence(tau_orient_seconds, orient_corr, ...)
     │
     ├── analysis.psd.ensemble_psd_analysis(data, dt)
     │    └── analysis.psd.compute_psd_components(r, dt)
     │         └── analysis.vacf.compute_velocity(r, dt)   [imported into psd]
     │
     ├── analysis.psd.fit_psd_powerlaw(freq_psd, ensemble_psd_vx/vy)
     │    └── analysis.psd._loglog_ols(logf, logp)
     │
     ├── analysis.diagnostics.retained_vs_discarded_lengths(lengths)
     ├── analysis.diagnostics.retained_vs_discarded_speed_stats(positions_all, lengths, dt)
     │    └── analysis.vacf.compute_velocity(r, dt)
     │
     └── analysis.plotting.*                               [all plot functions]
          (receive pre-computed arrays; no physics here)
```

---

## Shared dependency: `compute_velocity`

`analysis.vacf.compute_velocity` is the single canonical velocity function. It is imported by:

```
analysis.vacf          ← defines it
analysis.orientation   ← imports from analysis.vacf
analysis.psd           ← imports from analysis.vacf
analysis.diagnostics   ← imports from analysis.vacf
tracking.utils         ← imports from analysis.vacf (for CSV column only)
analysis.compute       ← re-exports as cal_velocity (shim)
```

---

## Shared dependency: `fit_msd_power_law`

Defined once in `analysis.msd`. Every code path that needs an MSD exponent uses this function:

```
analysis.msd.ensemble_msd_analysis()
    └── analysis.msd.fit_msd_power_law()   ← batch path

post_analysis.analysis_utils.fit_msd_power_law
    └── re-export of analysis.msd.fit_msd_power_law   ← notebook path

analysis.compute.fit_msd_power_law
    └── re-export of analysis.msd.fit_msd_power_law   ← legacy import path
```

All three resolve to the same function object. `alpha` is therefore bit-identical regardless of call site.

---

## `_loglog_ols` — OLS on log-log data

Two private copies exist (one in `msd.py`, one in `psd.py`) to avoid a circular import. Both are identical implementations. They:

1. Try `statsmodels.api.OLS` (full standard errors, R²)
2. Fall back to pure-NumPy OLS (same formulas) if `statsmodels` is absent

---

## Import dependency order (bottom → top)

```
numpy / scipy / matplotlib / statsmodels   (external)
    ↑
analysis.types
    ↑
analysis.vacf                  (compute_velocity defined here)
    ↑
analysis.msd                   (imports nothing from other analysis modules)
analysis.psd                   (imports compute_velocity from vacf)
analysis.orientation           (imports compute_velocity from vacf)
    ↑
analysis.diagnostics           (imports from msd and vacf)
    ↑
analysis.trajectories          (imports types only)
    ↑
analysis.plotting              (imports types only)
    ↑
analysis.main                  (imports from all of the above)
    ↑
analysis.compute               (shim — re-exports from above)
analysis.extract               (shim — re-exports from trajectories)
    ↑
tracking.utils                 (imports compute_velocity from vacf)
tracking.main                  (imports from tracking.utils, analysis.main)
    ↑
post_analysis.analysis_utils   (shim — re-exports from analysis.*)
    ↑
simulation.utils               (no analysis imports)
simulation.main                (imports simulation.utils)
    ↑
main.py                        (imports from simulation, tracking, analysis)
```

No circular imports exist.
