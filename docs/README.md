# Active-Droplet Project — Documentation

## What this project does

Simulates and tracks 2-D self-propelled active droplets, then applies a unified statistical analysis pipeline to characterise their dynamics.

Three independent entry paths exist:

| Path | What it does |
|------|-------------|
| `simulation/` | FEniCSx PDE + SDE simulation → trajectory CSV |
| `tracking/` | YOLO object detection on real videos → trajectory CSV |
| `analysis/` | Loads any trajectory CSV → MSD, VACF, PSD, orientation, plots |

All three paths produce or consume the **same CSV format** and feed into the **same canonical analysis functions**.

---

## Quick start

```python
from main import DropletPipeline

p = DropletPipeline()

# Run only analysis on existing CSVs
p.run_analysis(base_dir="output_droplet_tracking/csv_files",
               output_dir="./output",
               dt=0.2)

# Run video tracking + analysis
p.run_tracking(output_dir="output_droplet_tracking")

# Run simulation + (optionally) analysis
p.run_simulation()
```

---

## Documentation index

| Document | Contents |
|----------|----------|
| [Project structure](architecture/project_structure.md) | Every file and module explained |
| [Dependency graph](architecture/dependency_graph.md) | Full call chains |
| [Data flow](architecture/data_flow.md) | End-to-end pipeline diagram |
| [MSD](analysis/msd.md) | Mean squared displacement theory + implementation |
| [VACF](analysis/vacf.md) | Velocity autocorrelation + velocity definition |
| [Orientation](analysis/orientation.md) | Orientation correlation theory + implementation |
| [PSD](analysis/psd.md) | Power spectral density theory + implementation |
| [Persistence](analysis/persistence.md) | Translational and rotational persistence times |
| [Diagnostics](analysis/diagnostics.md) | Ensemble bias diagnostics |
| [Plotting conventions](analysis/plotting_conventions.md) | Standardized visual styling and deterministic coloring |
| [Tracking pipeline](tracking/tracking_pipeline.md) | YOLO tracking workflow |
| [Simulation overview](simulation/simulation_overview.md) | PDE + SDE simulation |
| [Function reference](api/function_reference.md) | Every public function |
| [Adding new analysis](developer/adding_new_analysis.md) | How to extend without breaking consistency |
| [Coding standards](developer/coding_standards.md) | Conventions enforced in this codebase |
| [Troubleshooting](developer/troubleshooting.md) | Debugging guide |

---

## CSV format

Every trajectory CSV has the following columns:

```
droplet_id, x, y, vx, vy, theta
```

- `droplet_id` — integer identifier
- `x`, `y` — position in pixels (tracking) or dimensionless units (simulation)
- `vx`, `vy` — forward-difference velocity (written at CSV creation time; **not used by analysis**)
- `theta` — heading angle in radians (written at CSV creation time; **not used by analysis**)

> **Note:** Analysis always recomputes velocity from `(x, y)` using the canonical forward-difference `v[i] = (r[i+1] − r[i]) / dt`. The stored `vx/vy/theta` columns are informational only.

---

## Python version and dependencies

See `requirements.txt`. Key dependencies:

| Package | Role |
|---------|------|
| `numpy` | Array operations |
| `scipy` | Welch PSD, curve fitting, optional Savitzky-Golay smoothing |
| `matplotlib` | All plotting |
| `statsmodels` | OLS regression for power-law fits (optional; NumPy fallback used if absent) |
| `ultralytics` | YOLO tracking (tracking path only) |
| `dolfinx / mpi4py / petsc4py` | FEniCSx PDE solver (simulation path only) |
| `pyvista` | 3-D visualisation for simulation GIF (simulation path only) |
