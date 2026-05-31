# Data Flow

End-to-end pipeline diagram using actual filenames and functions.

---

## Tracking path (real video → analysis)

```
videos/*.mp4
    │
    │  tracking.utils.generate_csvs()
    │  ultralytics YOLO .track()
    │  bounding box centres extracted
    │
    ▼
output_droplet_tracking/
    csv_files/<video_stem>_droplets.csv
    │   columns: droplet_id, x, y, vx, vy, theta
    │
    │  analysis.trajectories.load_positions_with_lengths(csv_path)
    │  → positions_consistent  Dict[int, ndarray(N,2)]
    │  → positions_all         Dict[int, ndarray(N_i,2)]
    │  → lengths               Dict[int, int]
    │
    ▼
positions_consistent  ← used for all ensemble quantities
positions_all         ← used for diagnostics (pair-weighted, speed stats)
    │
    │  analysis.msd.ensemble_msd_analysis()
    │  analysis.vacf.ensemble_vacf_analysis()
    │  analysis.orientation.ensemble_orientation_corr_analysis()
    │  analysis.psd.ensemble_psd_analysis()
    │  analysis.diagnostics.*
    │
    ▼
computed arrays  (tau, ensemble_msd, alpha, vacf_norm, orient_corr, ...)
    │
    │  analysis.plotting.*
    │
    ▼
output_droplet_tracking/
    plots/
        <label>_trajectories.png
        <label>_loglog_msd.png
        <label>_alpha_vs_tau.png
        <label>_msd_pair_weighted.png
        <label>_vacf.png
        <label>_orientation_corr.png
        <label>_psd_vx.png
        <label>_psd_vy.png
        <label>_track_length_hist.png
        <label>_results.txt
```

---

## Simulation path (PDE+SDE → analysis)

```
simulation.main.run_simulation()
    │
    │  FEniCSx mesh + function space
    │  Ndrop droplets, clustered initial positions
    │
    │  for each time step:
    │      update_paper_parameters(t, R0, Rshrink, ...)
    │      problem.solve()  ← diffusion PDE
    │      fem.assemble_scalar(Fx_form), Fy_form  ← gradient force
    │      x += -drift * Fx * dt + noise * randn  ← SDE update
    │      trajectory_x[i].append(Xnew)
    │
    ▼
simulation data/
    droplet_<N>_tracking.csv      ← same format as tracking CSV
    data_droplet_<N>_trajectories.npz
    multi_droplet_active.gif
    final_field.png
    droplet_<N>_trajectories.png
    │
    │  (same analysis pipeline as tracking path above)
    │  analysis.main.run_post_processing("simulation data", ...)
    │
    ▼
output_sim/   (or wherever output_dir points)
    <label>_loglog_msd.png
    <label>_results.txt
    ...
```

---

## Analysis-only path (existing CSVs → plots)

```
any directory containing *.csv files
    │
    │  DropletPipeline.run_analysis(base_dir, output_dir, dt)
    │  → analysis.main.run_post_processing(base_dir, output_dir, dt, ...)
    │
    ▼
(same as the lower half of the tracking path above)
```

---

## Data transformation detail

### CSV → `PositionsDict`

```python
# analysis/trajectories.py :: load_positions()

with open(csv_path, 'r') as f:
    reader = csv.reader(f)
    next(reader)                     # skip header row
    for row in reader:
        droplet_id = int(row[0])     # column 0
        x          = float(row[1])   # column 1
        y          = float(row[2])   # column 2
        # columns 3,4,5 (vx,vy,theta) are ignored

# Only full-length tracks are kept:
max_len = max(v.shape[0] for v in positions_by_id.values())
positions_consistent = {k: v for k, v in ... if v.shape[0] == max_len}
```

### `PositionsDict` → velocity

```python
# analysis/vacf.py :: compute_velocity()
v = (r[1:] - r[:-1]) / dt      # shape (N-1, 2)
```

This is the **only** velocity definition used in analysis. The `vx/vy` columns in the CSV are not read.

### Lag convention

```python
# Frame units (integers) — used internally:
tau = np.arange(1, max_lag + 1)      # for MSD
tau = np.arange(min_len)             # for VACF, orientation

# Seconds — used for axes and fitting:
tau_seconds = tau * dt               # in analysis/main.py
```

### Time step (`dt`) flow

`dt` is a caller-supplied parameter. It is not stored in the CSV.

| Path | How dt is set |
|------|--------------|
| Simulation | `dt_val = 0.002` (hardcoded in `simulation/main.py`) |
| Tracking | `dt = 1.0 / fps` where `fps = 30.0` (in `tracking/main.py`) |
| `run_analysis()` | Caller passes `dt`; default `0.2` in `DropletPipeline.run_analysis` |

> **Important:** If `dt` is wrong, all time-scaled quantities (τ_p, τ_r, frequency axis) will be wrong. The MSD **exponent** α is unaffected because it is computed on log(frame_index) vs log(MSD) and dt cancels out in the log-log slope. The local exponent uses `tau_seconds` as the x-axis, but d log(MSD) / d log(τ_seconds) = d log(MSD) / d log(tau * dt) = d log(MSD) / d log(tau) since dt is a constant offset in log space.

---

## Output files summary

| File | Produced by |
|------|-------------|
| `<label>_results.txt` | `analysis/main.py` |
| `<label>_loglog_msd.png` | `analysis/plotting.plot_msd_loglog` |
| `<label>_alpha_vs_tau.png` | `analysis/plotting.plot_local_msd_exponent` |
| `<label>_msd_pair_weighted.png` | `analysis/plotting.plot_msd_comparison` |
| `<label>_vacf.png` | `analysis/plotting.plot_vacf` |
| `<label>_orientation_corr.png` | `analysis/plotting.plot_orientation_corr` |
| `<label>_psd_vx.png` | `analysis/plotting.plot_psd_component` |
| `<label>_psd_vy.png` | `analysis/plotting.plot_psd_component` |
| `<label>_track_length_hist.png` | `analysis/plotting.plot_track_length_hist` |
| `<label>_trajectories.png` | `analysis/plotting.plot_trajectories` |
