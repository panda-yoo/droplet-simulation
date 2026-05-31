# Adding New Analysis

This guide explains how to extend the analysis package with new quantities without breaking the consistency guarantees established in the refactor.

---

## Rules that must not be violated

1. **One implementation per quantity.** Never write a function in a notebook or `post_analysis/` that computes something already in `analysis/`.
2. **All velocity must come from `analysis.vacf.compute_velocity`.** Do not write `(r[1:] - r[:-1]) / dt` anywhere else.
3. **Time convention:** Internal computations use integer frame indices. Convert to seconds with `tau_seconds = tau * dt` before plotting or fitting.
4. **Ensemble strategy:** Mean over per-entity estimates, after truncating to the shortest. Match the pattern in `ensemble_msd_analysis`, `ensemble_vacf_analysis`, etc.
5. **`analysis/plotting.py` is for plotting only.** Physics computations belong in the quantity module, not in plotting functions.

---

## Adding a new analysis module

### Step 1 — Create `analysis/<quantity>.py`

Use this template:

```python
"""
analysis.<quantity>
===================
Canonical <quantity> computations.

Contents
--------
- compute_<quantity>           : single-entity <quantity>
- ensemble_<quantity>_analysis : ensemble-averaged <quantity>
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple
from analysis.types import PositionsDict
from analysis.vacf import compute_velocity   # import velocity if needed


def compute_<quantity>(
    r: NDArray[np.floating],
    dt: float = 1.0,
    # ... parameters
) -> NDArray[np.floating]:
    """
    Compute <quantity> for a single trajectory.

    Parameters
    ----------
    r : ndarray, shape (N, 2)
        Ordered (x, y) positions.
    dt : float
        Time step between consecutive frames.

    Returns
    -------
    result : ndarray
        <quantity> at lags/times ... (frame units).
    """
    v = compute_velocity(r, dt)   # if velocity is needed
    # ... implementation


def ensemble_<quantity>_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]:
    """
    Compute ensemble-averaged <quantity>.

    Parameters
    ----------
    positions_dict : PositionsDict
        Mapping droplet_id -> (N, 2) position array.
    dt : float
        Time step.

    Returns
    -------
    tau : ndarray
        Lag times in frame units.
    ensemble_result : ndarray
        Mean over per-droplet <quantity>.
    """
    all_results = []
    for r in positions_dict.values():
        all_results.append(compute_<quantity>(r, dt))

    min_len = min(len(x) for x in all_results)
    ensemble_result = np.mean([x[:min_len] for x in all_results], axis=0)
    tau = np.arange(min_len)

    return tau, ensemble_result
```

### Step 2 — Add to `analysis/__init__.py`

```python
# In analysis/__init__.py, add a comment entry:
#   analysis.<quantity>  — canonical <quantity>
```

### Step 3 — Import in `analysis/main.py`

```python
# In analysis/main.py, add to imports:
from analysis.<quantity> import ensemble_<quantity>_analysis

# Inside post_processing(), add the computation:
tau_q, result_q = ensemble_<quantity>_analysis(data, dt)
tau_q_seconds = tau_q * dt
```

### Step 4 — Add a plot function in `analysis/plotting.py`

```python
def plot_<quantity>(
    tau_seconds: NDArray[np.floating],
    result: NDArray[np.floating],
    label: str,
    output_dir: str,
    n_ensemble: int,
) -> None:
    """Plot <quantity>."""
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(tau_seconds, result)
    ax.set_xlabel(r"Lag time $\tau$ (s)")
    ax.set_ylabel(r"$<quantity>$")
    ax.set_title(f"{label} : <quantity>\nN = {n_ensemble} droplets")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{label}_<quantity>.png"), dpi=300)
    plt.close(fig)
```

### Step 5 — Call the plot in `analysis/main.py`

```python
from analysis.plotting import plot_<quantity>

# Inside post_processing():
plot_<quantity>(tau_q_seconds, result_q, label, output_dir, n_ensemble)
```

### Step 6 — Add to results file

```python
# Inside the with open(...) block in post_processing():
f.write(f"<quantity> value = {result_q.mean():.6g}\n")
```

---

## Adding a new plot to an existing quantity

1. Add a new function to `analysis/plotting.py` following the existing pattern.
2. Call it from `analysis/main.py::post_processing` with pre-computed arrays.
3. Do not recompute physics inside the plot function.

---

## Adding a new ensemble quantity

If the new quantity follows the same pattern as MSD/VACF (loop over trajectories, compute per-entity, mean over ensemble), follow the template above exactly.

If the quantity has a different ensemble strategy (e.g. ensemble-level fit rather than mean), document the exception clearly in the module docstring.

---

## Adding a new diagnostic

1. Add a function to `analysis/diagnostics.py`.
2. Call it in `analysis/main.py::post_processing` with `positions_all` and/or `lengths`.
3. Write the result to `<label>_results.txt`.
4. Optionally add a plot.

Do not modify the retained/ensemble computation path. Diagnostics are read-only observers.

---

## Using the new quantity in a notebook

After completing steps 1–2, import directly from the canonical module:

```python
from analysis.<quantity> import compute_<quantity>, ensemble_<quantity>_analysis
from analysis.trajectories import load_positions
from pathlib import Path

positions, _, _ = load_positions(Path("data.csv"))
dt = 0.2

tau, result = ensemble_<quantity>_analysis(positions, dt)
tau_seconds = tau * dt
```

Do not import from `post_analysis/` for new quantities.

---

## Checklist before committing

- [ ] New function has a NumPy-style docstring (Parameters / Returns / Notes)
- [ ] All parameters have type hints
- [ ] Velocity uses `from analysis.vacf import compute_velocity`
- [ ] Lag arrays returned in frame units; seconds conversion done by caller
- [ ] Ensemble uses `np.mean([x[:min_len] for x in all_results], axis=0)`
- [ ] Plot function only plots — no physics
- [ ] `analysis/compute.py` shim updated if backward compat is needed
