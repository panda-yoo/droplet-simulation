# Coding Standards

Standards enforced throughout this codebase. Follow these when adding or modifying any module.

---

## 1. One implementation per quantity

Every physical quantity has exactly one canonical implementation:

| Quantity | Canonical location |
|----------|-------------------|
| Velocity | `analysis/vacf.py::compute_velocity` |
| MSD | `analysis/msd.py::compute_msd` |
| MSD fit | `analysis/msd.py::fit_msd_power_law` |
| Local MSD exponent | `analysis/msd.py::compute_local_msd_exponent` |
| VACF | `analysis/vacf.py::compute_vacf` |
| Orientation angle | `analysis/orientation.py::compute_heading_angle` |
| Orientation correlation | `analysis/orientation.py::compute_orientation_corr` |
| PSD | `analysis/psd.py::compute_psd_components` |

Do not write an alternative implementation anywhere. If a notebook needs MSD, import `compute_msd` from `analysis.msd`.

---

## 2. Type hints on all public functions

Every public function must have:
- Type hints on all parameters
- Type hint on return value

```python
# Correct
def compute_msd(
    r: NDArray[np.floating],
    lag_fraction: float = 1.0,
    minimum_pairs: Optional[int] = None,
) -> NDArray[np.floating]:

# Wrong — missing type hints
def compute_msd(r, lag_fraction=1.0):
```

Use `NDArray[np.floating]` for float arrays, `NDArray[np.int_]` for integer arrays, `PositionsDict` for trajectory dicts.

---

## 3. NumPy-style docstrings on all public functions

Every public function must have a docstring with at minimum `Parameters` and `Returns` sections:

```python
def compute_msd(
    r: NDArray[np.floating],
    lag_fraction: float = 1.0,
) -> NDArray[np.floating]:
    """
    One-line summary.

    Parameters
    ----------
    r : ndarray, shape (N, 2)
        Description of r.
    lag_fraction : float
        Description.

    Returns
    -------
    msd : ndarray, shape (max_lag,)
        Description.

    Notes
    -----
    Optional: algorithm details, assumptions, references.
    """
```

Internal helper functions (prefixed with `_`) do not require full docstrings but should have a one-line description.

---

## 4. Time convention

**Frame units internally:**
```python
tau = np.arange(1, max_lag + 1)    # integers, frame units
```

**Seconds for plotting and fitting:**
```python
tau_seconds = tau * dt             # in analysis/main.py, not in the quantity module
```

Quantity modules return frame-unit lags. Callers convert. The quantity modules must not assume a value of `dt` in their output arrays.

---

## 5. Lag indexing convention

- **MSD:** `tau[i]` corresponds to lag `τ = i + 1` (because lag 0 is trivially 0)
- **VACF, orientation:** `tau[i]` corresponds to lag `τ = i` (lag 0 is the self-correlation)

This is already established. Do not change it.

---

## 6. Ensemble averaging pattern

All ensemble functions follow the same pattern:

```python
all_results = [compute_X(r, ...) for r in positions_dict.values()]
min_len = min(len(x) for x in all_results)
ensemble_result = np.mean([x[:min_len] for x in all_results], axis=0)
```

Do not use `np.mean` on ragged arrays or pad short arrays with zeros.

---

## 7. Plotting functions do not compute physics

Functions in `analysis/plotting.py` receive pre-computed arrays. They do not call any function from `msd.py`, `vacf.py`, etc. This ensures that the number plotted equals the number computed, with no possibility of divergence.

```python
# Correct
def plot_msd_loglog(tau_seconds, ensemble_msd, alpha, fit_line_full, ...):
    ax.loglog(tau_seconds, ensemble_msd)
    # alpha was computed elsewhere and passed in

# Wrong — physics inside plotting
def plot_msd_loglog(positions_dict, dt, ...):
    _, msd, alpha, *_ = ensemble_msd_analysis(positions_dict)  # BAD
```

---

## 8. `statsmodels` optional dependency

`statsmodels` provides standard errors for OLS regression. When it is not installed, a pure-NumPy fallback is used. Both produce identical point estimates (slope = alpha or beta); only the standard errors may differ slightly.

Pattern:

```python
try:
    import statsmodels.api as sm
except Exception:
    sm = None

def _loglog_ols(...):
    if sm is not None:
        # use statsmodels
    else:
        # pure-NumPy fallback
```

Do not make `statsmodels` a hard requirement in any new code.

---

## 9. `scipy.optimize.curve_fit` must be guarded

```python
try:
    popt, pcov = curve_fit(...)
except RuntimeError:
    return np.nan, np.nan, tau_fit, np.array([])
```

Fitting can fail if the initial guess is poor or data is noisy. Never let a fit failure crash the pipeline. Return NaN values on failure.

---

## 10. Compatibility shims

`analysis/compute.py`, `analysis/extract.py`, and `post_analysis/analysis_utils.py` are compatibility shims. They re-export canonical names so that existing notebooks continue to work.

Rules:
- **Do not add new implementations to shims.**
- When you add a new function to a canonical module and you expect notebooks to import it from the shim, add the re-export.
- Document shim additions with a comment: `# Added for backward compat`.

---

## 11. File naming

| Content | Module |
|---------|--------|
| CSV loading | `analysis/trajectories.py` |
| MSD | `analysis/msd.py` |
| VACF + velocity | `analysis/vacf.py` |
| Orientation | `analysis/orientation.py` |
| PSD | `analysis/psd.py` |
| Diagnostics | `analysis/diagnostics.py` |
| Plotting | `analysis/plotting.py` |
| Orchestration | `analysis/main.py` |
| Type aliases | `analysis/types.py` |

New quantities get their own file, not a section in `compute.py`.

---

## 12. No global mutable state

Do not use module-level mutable variables. Every function receives all parameters it needs via arguments. This ensures that calling the same function twice with the same arguments always returns the same result.
