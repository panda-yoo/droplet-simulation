# Project Structure

## Repository root

```
combine - Copy/
├── main.py                     ← Unified pipeline entry point (DropletPipeline class)
├── requirements.txt
├── .gitignore
│
├── analysis/                   ← Canonical analysis package
├── post_analysis/              ← Legacy namespace (shims only — no implementations)
├── tracking/                   ← YOLO tracking pipeline
├── simulation/                 ← FEniCSx PDE + SDE simulation
│
├── docs/                       ← This documentation
├── videos/                     ← Input videos for tracking
├── weight/                     ← YOLO model weights (best.pt)
├── output/                     ← Default analysis output directory
├── output_droplet_tracking/    ← Default tracking output directory
│
├── nb1.ipynb                   ← Notebook (analysis / plotting)
└── nb2.ipynb                   ← Notebook (analysis / plotting)
```

---

## `main.py` — Unified entry point

**Purpose:** Single class (`DropletPipeline`) that exposes `run_simulation()`, `run_tracking()`, `run_analysis()`, and `run_all()`.

**Responsibilities:**
- Route calls to the correct sub-package
- Expose a clean API for notebook and CLI use
- No physics or analysis logic lives here

**Public class:** `DropletPipeline`

**Called by:** notebooks, CLI (`python main.py`)

**Calls:**
- `simulation.main.run_simulation`
- `tracking.main.run_tracking`
- `analysis.main.run_post_processing`

---

## `analysis/` — Canonical analysis package

All physics computations live here. Every public function is documented, type-hinted, and has a NumPy-style docstring.

### `analysis/types.py`

**Purpose:** Shared type aliases.

**Public names:**
- `PositionsDict = Dict[int, NDArray[np.floating]]` — the canonical type for a collection of trajectories

**Used by:** every analysis module

---

### `analysis/trajectories.py`

**Purpose:** Load droplet trajectory data from CSV files.

**Responsibilities:**
- Parse CSV rows into per-droplet numpy arrays
- Apply the full-length filter (discard partial tracks)
- Return consistent and all-track dicts

**Public functions:**

| Function | Returns |
|----------|---------|
| `load_positions(csv_path)` | `(positions_consistent, total_paths, consistent_paths)` |
| `load_positions_with_lengths(csv_path)` | above + `positions_all`, `lengths` |

**Dependencies:** `csv`, `numpy`, `analysis.types`

**Used by:** `analysis/main.py`, `post_analysis/analysis_utils.py` (shim), notebooks

---

### `analysis/msd.py`

**Purpose:** Single canonical MSD implementation.

**Responsibilities:**
- Single-droplet MSD computation
- Power-law fitting on log-log scale
- Ensemble-averaged MSD + exponent (delegates to `fit_msd_power_law`)
- Local (scale-dependent) MSD exponent
- Pair-weighted MSD for bias diagnosis

**Public functions:**

| Function | Description |
|----------|-------------|
| `compute_msd(r, lag_fraction, minimum_pairs)` | Single-droplet MSD |
| `fit_msd_power_law(msd, fit_fraction, minimum_pairs)` | OLS fit MSD ~ τ^α |
| `ensemble_msd_analysis(positions_dict, ...)` | Ensemble MSD + α |
| `compute_local_msd_exponent(tau, msd)` | d log MSD / d log τ |
| `compute_pair_weighted_msd(positions_dict, ...)` | Pair-weighted MSD |

**Internal:** `_loglog_ols` — OLS on log-log data (uses `statsmodels` when available)

**Dependencies:** `numpy`, `statsmodels` (optional), `analysis.types`

**Used by:** `analysis/main.py`, `analysis/plotting.py`, `analysis/diagnostics.py`, `post_analysis/analysis_utils.py` (shim), notebooks

---

### `analysis/vacf.py`

**Purpose:** Single canonical VACF implementation. Also defines the **sole velocity function** used by the entire analysis package.

**Responsibilities:**
- Forward-difference velocity (the canonical definition)
- Single-droplet VACF
- Exponential persistence-time fit
- Ensemble-averaged VACF

**Public functions:**

| Function | Description |
|----------|-------------|
| `compute_velocity(r, dt)` | `v = (r[1:] - r[:-1]) / dt`, shape `(N-1, 2)` |
| `compute_vacf(r, dt)` | Single-droplet `C_v(τ)` |
| `fit_vacf_exponential(tau, vacf, max_fit_fraction)` | Fit `A exp(-t/τ_p)` |
| `ensemble_vacf_analysis(positions_dict, dt)` | Ensemble VACF raw + normalised |

**Dependencies:** `numpy`, `scipy.optimize.curve_fit`, `analysis.types`

**Used by:** `analysis/main.py`, `analysis/orientation.py` (imports `compute_velocity`), `analysis/psd.py` (imports `compute_velocity`), `analysis/diagnostics.py` (imports `compute_velocity`), `tracking/utils.py` (imports `compute_velocity`)

---

### `analysis/orientation.py`

**Purpose:** Single canonical orientation correlation implementation.

**Responsibilities:**
- Heading angle from velocity (forward difference, unwrapped)
- Single-droplet `⟨cos(Δθ(τ))⟩` correlation
- Rotational persistence-time fit
- Ensemble-averaged orientation correlation

**Public functions:**

| Function | Description |
|----------|-------------|
| `compute_heading_angle(r, dt)` | Unwrapped `arctan2(vy, vx)`, shape `(N-1,)` |
| `compute_orientation_corr(theta)` | `C(τ) = ⟨cos(θ(t+τ) − θ(t))⟩` |
| `fit_orientation_persistence(tau, corr, max_fit_fraction)` | Fit `exp(-t/τ_r)` |
| `ensemble_orientation_corr_analysis(positions_dict, dt)` | Ensemble orientation corr |

**Dependencies:** `numpy`, `scipy.optimize.curve_fit`, `analysis.vacf.compute_velocity`, `analysis.types`

**Used by:** `analysis/main.py`, `analysis/plotting.py`

---

### `analysis/psd.py`

**Purpose:** Single canonical PSD implementation using Welch's method.

**Responsibilities:**
- Welch PSD of vx and vy for a single trajectory
- Power-law fit PSD ~ f^{−β}
- Ensemble-averaged PSD

**Public functions:**

| Function | Description |
|----------|-------------|
| `compute_psd_components(r, dt)` | Welch PSD of mean-subtracted vx, vy |
| `fit_psd_powerlaw(freqs, psd, fmin, fmax)` | OLS fit PSD ~ f^{−β} |
| `ensemble_psd_analysis(positions_dict, dt)` | Ensemble-averaged PSD |

**Internal:** `_loglog_ols` — same OLS routine as in `msd.py` (local copy to avoid circular import)

**Dependencies:** `numpy`, `scipy.signal.welch`, `statsmodels` (optional), `analysis.vacf.compute_velocity`, `analysis.types`

**Used by:** `analysis/main.py`, `analysis/plotting.py`

---

### `analysis/diagnostics.py`

**Purpose:** Ensemble bias diagnostics. Quantify the effect of the full-length track filter.

**Responsibilities:**
- Track-length distribution (retained vs discarded)
- Speed statistics (retained vs discarded)
- Pair-weighted MSD (re-export from `msd.py`)
- Pair-weighted VACF

**Public functions:**

| Function | Description |
|----------|-------------|
| `retained_vs_discarded_lengths(lengths)` | Split lengths by retained/discarded |
| `retained_vs_discarded_speed_stats(positions_all, lengths, dt)` | Speed stats |
| `pair_weighted_msd` | Re-export of `compute_pair_weighted_msd` |
| `pair_weighted_vacf(positions_dict, dt, lag_fraction)` | Pair-weighted VACF |

**Dependencies:** `numpy`, `analysis.msd`, `analysis.vacf`, `analysis.types`

**Used by:** `analysis/main.py`

---

### `analysis/plotting.py`

**Purpose:** All plot generation. No physics is computed here.

**Responsibilities:**
- Accept pre-computed quantities and save PNG figures
- Enforce consistent axis labels and time conventions

**Public functions:**

| Function | Output file |
|----------|-------------|
| `plot_trajectories(...)` | `<label>_trajectories.png` |
| `plot_msd_loglog(...)` | `<label>_loglog_msd.png` |
| `plot_vacf(...)` | `<label>_vacf.png` |
| `plot_orientation_corr(...)` | `<label>_orientation_corr.png` |
| `plot_psd_component(...)` | `<label>_psd_vx.png` / `<label>_psd_vy.png` |
| `plot_local_msd_exponent(...)` | `<label>_alpha_vs_tau.png` |
| `plot_msd_comparison(...)` | `<label>_msd_pair_weighted.png` |
| `plot_track_length_hist(...)` | `<label>_track_length_hist.png` |

**Dependencies:** `numpy`, `matplotlib`, `analysis.types`

**Used by:** `analysis/main.py`

---

### `analysis/main.py`

**Purpose:** Orchestration. Loads data, calls canonical analysis functions, calls plotting functions, writes results files.

**Responsibilities:**
- Iterate over CSVs in a folder
- Call every canonical ensemble function
- Convert frame-index lags to seconds (`tau_seconds = tau * dt`)
- Delegate all plotting to `analysis/plotting.py`
- Write `<label>_results.txt`

**Public functions:**

| Function | Description |
|----------|-------------|
| `post_processing(folder_path, output_dir, dt, ...)` | Process all CSVs in one folder |
| `run_post_processing(base_dir, output_dir, dt, ...)` | Find CSVs and call `post_processing` |

**Dependencies:** All analysis modules, `analysis.trajectories`, `os`, `numpy`

**Called by:** `main.py`, `tracking/main.py`

---

### `analysis/compute.py` — Compatibility shim

**Purpose:** Re-exports every name that was previously defined in this file.  
**Do not add implementations here.**

Maps: `cal_velocity` → `compute_velocity`, all MSD/VACF/PSD/orientation functions → canonical modules.

---

### `analysis/extract.py` — Compatibility shim

**Purpose:** Re-exports `load_positions` and `load_positions_with_lengths` from `analysis/trajectories.py`.

---

### `analysis/__init__.py`

**Purpose:** Package init. Imports `load_positions` and `load_positions_with_lengths` for convenience.

---

## `post_analysis/` — Legacy namespace (shims only)

### `post_analysis/analysis_utils.py` — Compatibility shim

**Purpose:** Re-exports every public name from the canonical `analysis.*` modules.  
Notebooks that imported `from post_analysis.analysis_utils import fit_msd_power_law` continue to work and now call the same canonical function as the batch pipeline.

**Do not add implementations here.**

---

## `tracking/`

### `tracking/main.py`

**Purpose:** YOLO tracking entry point.

**Responsibilities:**
- Load YOLO model from `weight/best.pt`
- Discover video files in `videos/`
- Call `generate_csvs` to produce per-video CSVs
- Call `analysis.main.run_post_processing` to analyse them

**Public function:** `run_tracking(output_dir, save_annotated)`

**Dependencies:** `ultralytics`, `tracking.utils`, `analysis.main`

---

### `tracking/utils.py`

**Purpose:** YOLO tracking helpers.

**Responsibilities:**
- Convert YOLO bounding boxes to centre-point positions
- Compute velocity for CSV columns only (padded forward difference via `analysis.vacf.compute_velocity`)
- Write per-video trajectory CSVs

**Public functions:**

| Function | Description |
|----------|-------------|
| `positions_to_rows(positions_by_id, dt)` | Convert dict → flat CSV rows |
| `generate_csvs(model, videos_dir, ...)` | Run YOLO, collect positions, write CSVs |

**Internal:** `_velocity_for_csv(pos, dt)` — padded velocity for CSV column only; uses `compute_velocity` internally

**Dependencies:** `ultralytics`, `numpy`, `csv`, `analysis.vacf.compute_velocity`

---

## `simulation/`

### `simulation/main.py`

**Purpose:** Full 2-D active-droplet simulation.

**Physics:** Chemical trail diffusion PDE (FEniCSx) + stochastic SDE for droplet motion.

**Responsibilities:**
- Build FE mesh and function space (DOLFINx)
- Place droplets at clustered initial positions
- Time-step: solve diffusion PDE → assemble gradient force → SDE update
- Record animated GIF (PyVista)
- Save trajectories (.npz and .csv)

**Public function:** `run_simulation()`

**Dependencies:** `dolfinx`, `mpi4py`, `petsc4py`, `pyvista`, `numpy`, `matplotlib`, `tqdm`, `simulation.utils`

---

### `simulation/utils.py`

**Purpose:** Helper functions for the simulation.

**Public functions:**

| Function | Description |
|----------|-------------|
| `initial_condition(x)` | Zero-field initial condition for the chemical PDE |
| `update_paper_parameters(t, R0, Rshrink, Rmin, delta, eta)` | Time-dependent R, α, γ |
| `save_trajectory_csv(trajectory_x, trajectory_y, dt_val, Ndrop, output_path)` | Write trajectory CSV |

**Dependencies:** `numpy`, `csv`
