# Function Reference

Complete reference for every public function in the analysis package, generated from the actual source code.

---

## `analysis.types`

### `PositionsDict`

```python
PositionsDict = Dict[int, NDArray[np.floating]]
```

Type alias for the canonical trajectory collection: maps integer droplet ID to a float array of shape `(N, 2)` containing ordered `(x, y)` positions.

---

## `analysis.trajectories`

### `load_positions`

```python
def load_positions(csv_path: Path) -> Tuple[PositionsDict, int, int]
```

| | |
|--|--|
| **Location** | `analysis/trajectories.py` |
| **Parameters** | `csv_path: Path` — path to tracking CSV |
| **Returns** | `(positions_consistent, total_paths, consistent_paths)` |
| `positions_consistent` | Full-length tracks only, `Dict[int, ndarray(N,2)]` |
| `total_paths` | Total unique droplet IDs |
| `consistent_paths` | Number of retained IDs |
| **Called by** | `load_positions_with_lengths`, `post_analysis/analysis_utils.individual_analysis` |
| **Notes** | Only keeps tracks with length equal to the global maximum. Ignores `vx`, `vy`, `theta` CSV columns. |

---

### `load_positions_with_lengths`

```python
def load_positions_with_lengths(csv_path: Path) -> Tuple[
    PositionsDict, int, int, PositionsDict, Dict[int, int]
]
```

| | |
|--|--|
| **Location** | `analysis/trajectories.py` |
| **Parameters** | `csv_path: Path` |
| **Returns** | `(positions_consistent, total_paths, consistent_paths, positions_all, lengths)` |
| `positions_all` | All trajectories including partial tracks |
| `lengths` | `Dict[int, int]` — frames per droplet ID |
| **Called by** | `analysis/main.py::post_processing` |

---

## `analysis.msd`

### `compute_msd`

```python
def compute_msd(
    r: NDArray[np.floating],
    lag_fraction: float = 1.0,
    minimum_pairs: Optional[int] = None,
) -> NDArray[np.floating]
```

| | |
|--|--|
| **Location** | `analysis/msd.py` |
| **Parameters** | `r` — `(N, 2)` positions; `lag_fraction` — max lag as fraction of N-1; `minimum_pairs` — minimum displacement pairs at any lag |
| **Returns** | `ndarray(max_lag,)` — MSD at lags τ=1,…,max_lag (frame units) |
| **Called by** | `ensemble_msd_analysis`, `post_analysis/analysis_utils` (shim) |
| **Complexity** | O(N × max_lag) |
| **Notes** | Direct sliding-window estimator. Index `i` corresponds to lag τ=i+1. |

---

### `fit_msd_power_law`

```python
def fit_msd_power_law(
    msd: NDArray[np.floating],
    fit_fraction: float = 0.3,
    minimum_pairs: Optional[int] = None,
) -> Tuple[float, float, float, NDArray[np.floating]]
```

| | |
|--|--|
| **Location** | `analysis/msd.py` |
| **Parameters** | `msd` — MSD array (frame units); `fit_fraction` — fraction of MSD length to fit; `minimum_pairs` — additional lag cap |
| **Returns** | `(alpha, alpha_err, r2, fit_line_full)` |
| `alpha` | Power-law exponent (slope of log10 MSD vs log10 τ) |
| `alpha_err` | Standard error of alpha |
| `r2` | Coefficient of determination |
| `fit_line_full` | Fitted values aligned with full MSD array; NaN outside fit region |
| **Called by** | `ensemble_msd_analysis`, `post_analysis/analysis_utils` (shim) |
| **Notes** | Uses `statsmodels.OLS` when available, pure-NumPy OLS fallback otherwise. Both produce identical results. |

---

### `ensemble_msd_analysis`

```python
def ensemble_msd_analysis(
    positions_dict: PositionsDict,
    lag_fraction: float = 0.3,
    fit_fraction: float = 0.3,
    minimum_pairs: Optional[int] = None,
) -> Tuple[NDArray[np.int_], NDArray[np.floating], float, float, float, NDArray[np.floating]]
```

| | |
|--|--|
| **Location** | `analysis/msd.py` |
| **Parameters** | `positions_dict` — trajectory collection; `lag_fraction` — max lag fraction; `fit_fraction` — fit region fraction; `minimum_pairs` — pair count floor |
| **Returns** | `(tau, ensemble_msd, alpha, alpha_err, r2, fit_line_full)` |
| `tau` | Lag times in frame units: 1, 2, …, L |
| `ensemble_msd` | Mean over per-droplet MSDs |
| **Called by** | `analysis/main.py::post_processing` |
| **Notes** | Internally calls `compute_msd` then `fit_msd_power_law`. Alpha is therefore identical to what `fit_msd_power_law` would return on the same ensemble MSD. |

---

### `compute_local_msd_exponent`

```python
def compute_local_msd_exponent(
    tau: NDArray[np.floating],
    msd: NDArray[np.floating],
) -> NDArray[np.floating]
```

| | |
|--|--|
| **Location** | `analysis/msd.py` |
| **Parameters** | `tau` — positive strictly-increasing lags (any units); `msd` — positive MSD values |
| **Returns** | `ndarray` — d log10(MSD) / d log10(τ) at every lag, via central differences |
| **Called by** | `analysis/main.py::post_processing` |
| **Notes** | Raises `ValueError` if inputs are non-finite, non-positive, or non-monotone. Unit-agnostic (pass frame integers or seconds — the log ratio is the same). |

---

### `compute_pair_weighted_msd`

```python
def compute_pair_weighted_msd(
    positions_dict: PositionsDict,
    lag_fraction: float = 0.3,
    minimum_pairs: Optional[int] = None,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]
```

| | |
|--|--|
| **Location** | `analysis/msd.py` |
| **Parameters** | `positions_dict` — trajectories (can include partial tracks); `lag_fraction` — based on shortest track |
| **Returns** | `(tau, msd_weighted)` — frame-unit lags and pair-weighted MSD |
| **Called by** | `analysis/main.py`, `analysis/diagnostics.pair_weighted_msd` (re-export) |
| **Notes** | Each displacement pair contributes equally regardless of which trajectory it came from. |

---

## `analysis.vacf`

### `compute_velocity`

```python
def compute_velocity(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]
```

| | |
|--|--|
| **Location** | `analysis/vacf.py` |
| **Parameters** | `r` — `(N, 2)` positions; `dt` — time step |
| **Returns** | `ndarray(N-1, 2)` — forward-difference velocity |
| **Called by** | `compute_vacf`, `orientation.compute_heading_angle`, `psd.compute_psd_components`, `diagnostics.retained_vs_discarded_speed_stats`, `diagnostics.pair_weighted_vacf`, `tracking/utils._velocity_for_csv` |
| **Notes** | `v[i] = (r[i+1] - r[i]) / dt`. This is the **only** velocity definition used in the project. |

---

### `compute_vacf`

```python
def compute_vacf(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]
```

| | |
|--|--|
| **Location** | `analysis/vacf.py` |
| **Parameters** | `r` — `(N, 2)` positions; `dt` — time step |
| **Returns** | `ndarray(N-1,)` — VACF at lags 0, 1, …, N-2 (frame units) |
| **Called by** | `ensemble_vacf_analysis` |
| **Complexity** | O((N-1)²) |

---

### `fit_vacf_exponential`

```python
def fit_vacf_exponential(
    tau: NDArray[np.floating],
    vacf: NDArray[np.floating],
    max_fit_fraction: Optional[float] = None,
) -> Tuple[float, float, NDArray[np.floating], NDArray[np.floating]]
```

| | |
|--|--|
| **Location** | `analysis/vacf.py` |
| **Parameters** | `tau` — lag times (any units); `vacf` — VACF values; `max_fit_fraction` — optional cap |
| **Returns** | `(tau_p, tau_p_err, tau_fit, fit_curve)` |
| `tau_p` | Velocity persistence time (same units as `tau`) |
| `tau_p_err` | Standard error of τ_p |
| **Called by** | `analysis/main.py::post_processing` |
| **Notes** | Stops fit at first zero crossing. Returns NaN if fewer than 5 fit points. Uses `scipy.optimize.curve_fit`. |

---

### `ensemble_vacf_analysis`

```python
def ensemble_vacf_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray[np.int_], NDArray[np.floating], NDArray[np.floating]]
```

| | |
|--|--|
| **Location** | `analysis/vacf.py` |
| **Parameters** | `positions_dict` — trajectories; `dt` — time step |
| **Returns** | `(tau, vacf_raw, vacf_norm)` — frame-unit lags, raw VACF, normalised VACF |
| **Called by** | `analysis/main.py::post_processing` |

---

## `analysis.orientation`

### `compute_heading_angle`

```python
def compute_heading_angle(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]
```

| | |
|--|--|
| **Location** | `analysis/orientation.py` |
| **Parameters** | `r` — `(N, 2)` positions; `dt` — time step |
| **Returns** | `ndarray(N-1,)` — unwrapped heading angle (radians) |
| **Called by** | `ensemble_orientation_corr_analysis` |

---

### `compute_orientation_corr`

```python
def compute_orientation_corr(
    theta: NDArray[np.floating],
) -> NDArray[np.floating]
```

| | |
|--|--|
| **Location** | `analysis/orientation.py` |
| **Parameters** | `theta` — heading angles (radians) |
| **Returns** | `ndarray(M,)` — `⟨cos(θ(t+τ)−θ(t))⟩` at lags 0,…,M-1 |
| **Called by** | `ensemble_orientation_corr_analysis` |
| **Complexity** | O(M²) |

---

### `fit_orientation_persistence`

```python
def fit_orientation_persistence(
    tau: NDArray[np.floating],
    corr: NDArray[np.floating],
    max_fit_fraction: Optional[float] = None,
) -> Tuple[float, float, NDArray[np.floating], NDArray[np.floating]]
```

| | |
|--|--|
| **Location** | `analysis/orientation.py` |
| **Parameters** | `tau` — lags; `corr` — orientation correlation; `max_fit_fraction` — optional cap |
| **Returns** | `(tau_r, tau_r_err, tau_fit, fit_curve)` |
| `tau_r` | Rotational persistence time (same units as `tau`) |
| **Called by** | `analysis/main.py::post_processing` |
| **Notes** | Model: `exp(-t/τ_r)` (no amplitude — C_θ(0)=1). Stops at first zero crossing. |

---

### `ensemble_orientation_corr_analysis`

```python
def ensemble_orientation_corr_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]
```

| | |
|--|--|
| **Location** | `analysis/orientation.py` |
| **Parameters** | `positions_dict` — trajectories; `dt` — time step |
| **Returns** | `(tau, ensemble_corr)` — frame-unit lags, ensemble `C_θ(τ)` |
| **Called by** | `analysis/main.py::post_processing` |

---

## `analysis.psd`

### `compute_psd_components`

```python
def compute_psd_components(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> Tuple[NDArray, NDArray, NDArray]
```

| | |
|--|--|
| **Location** | `analysis/psd.py` |
| **Parameters** | `r` — `(N, 2)` positions; `dt` — time step |
| **Returns** | `(freq, psd_vx, psd_vy)` — frequencies in 1/dt units, Welch PSDs |
| **Called by** | `ensemble_psd_analysis` |
| **Notes** | Returns three empty arrays if `len(v) < 8`. `nperseg = min(256, len(v))`. |

---

### `fit_psd_powerlaw`

```python
def fit_psd_powerlaw(
    freqs: NDArray[np.floating],
    psd: NDArray[np.floating],
    fmin: float = 0.1,
    fmax: float = 2.0,
) -> Tuple[float, float, float, NDArray[np.floating], NDArray[np.floating]]
```

| | |
|--|--|
| **Location** | `analysis/psd.py` |
| **Parameters** | `freqs` — frequency array; `psd` — PSD values; `fmin`, `fmax` — fit frequency bounds |
| **Returns** | `(beta, beta_err, r2, fit_freqs, fit_line)` |
| `beta` | Positive exponent: PSD ~ f^{-β} |
| **Called by** | `analysis/main.py::post_processing` |
| **Notes** | Returns NaN if fewer than 3 points in [fmin, fmax]. |

---

### `ensemble_psd_analysis`

```python
def ensemble_psd_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray, NDArray, NDArray]
```

| | |
|--|--|
| **Location** | `analysis/psd.py` |
| **Parameters** | `positions_dict` — trajectories; `dt` — time step |
| **Returns** | `(freq, ensemble_psd_vx, ensemble_psd_vy)` |
| **Called by** | `analysis/main.py::post_processing` |

---

## `analysis.diagnostics`

### `retained_vs_discarded_lengths`

```python
def retained_vs_discarded_lengths(
    lengths: Dict[int, int],
) -> Tuple[List[int], List[int], int, int]
```

| | |
|--|--|
| **Location** | `analysis/diagnostics.py` |
| **Parameters** | `lengths` — track lengths from `load_positions_with_lengths` |
| **Returns** | `(retained_lengths, discarded_lengths, n_retained, n_discarded)` |
| **Called by** | `analysis/main.py::post_processing` |

---

### `retained_vs_discarded_speed_stats`

```python
def retained_vs_discarded_speed_stats(
    positions_all: PositionsDict,
    lengths: Dict[int, int],
    dt: float = 1.0,
) -> Dict[str, float]
```

| | |
|--|--|
| **Location** | `analysis/diagnostics.py` |
| **Parameters** | `positions_all` — all tracks; `lengths` — track lengths; `dt` — time step |
| **Returns** | Dict with `retained_mean_speed`, `retained_std_speed`, `discarded_mean_speed`, `discarded_std_speed`, `n_retained`, `n_discarded` |
| **Called by** | `analysis/main.py::post_processing` |

---

### `pair_weighted_msd`

Re-export of `analysis.msd.compute_pair_weighted_msd`. See above.

---

### `pair_weighted_vacf`

```python
def pair_weighted_vacf(
    positions_dict: PositionsDict,
    dt: float = 1.0,
    lag_fraction: float = 0.3,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]
```

| | |
|--|--|
| **Location** | `analysis/diagnostics.py` |
| **Parameters** | `positions_dict` — any trajectory collection; `dt` — time step; `lag_fraction` — fraction of shortest track |
| **Returns** | `(tau, vacf_weighted)` — frame-unit lags, pair-weighted VACF |
| **Called by** | Not called by default pipeline; import directly for manual bias analysis |

---

## `analysis.main`

### `post_processing`

```python
def post_processing(
    folder_path: str,
    output_dir: str,
    dt: float = 1.0,
    msd_lag_fraction: float = 0.3,
    msd_fit_fraction: float = 0.3,
    msd_minimum_pairs: Optional[int] = None,
    vacf_max_fit_fraction: Optional[float] = None,
    orient_max_fit_fraction: Optional[float] = None,
) -> None
```

| | |
|--|--|
| **Location** | `analysis/main.py` |
| **Parameters** | See docstring |
| **Returns** | None (saves files to `output_dir`) |
| **Called by** | `run_post_processing`, `tracking/main.py` |

---

### `run_post_processing`

```python
def run_post_processing(
    base_dir: str,
    output_dir: str,
    dt: float = 1.0,
    ...
) -> None
```

| | |
|--|--|
| **Location** | `analysis/main.py` |
| **Parameters** | Same as `post_processing` plus `base_dir` |
| **Called by** | `main.py::DropletPipeline.run_analysis`, `tracking/main.py::run_tracking` |
| **Notes** | Processes CSVs directly in `base_dir` and in all sub-folders. |

---

## `simulation.utils`

### `update_paper_parameters`

```python
def update_paper_parameters(t, R0, Rshrink, Rmin, delta, eta) -> Tuple[float, float, float]
```

| | |
|--|--|
| **Location** | `simulation/utils.py` |
| **Returns** | `(Rval, alphaval, gammaval)` — radius, secretion rate, Stokes friction |
| **Called by** | `simulation/main.py` at every time step |

---

### `save_trajectory_csv`

```python
def save_trajectory_csv(
    trajectory_x: List[List[float]],
    trajectory_y: List[List[float]],
    dt_val: float,
    Ndrop: int,
    output_path: str | Path,
) -> None
```

| | |
|--|--|
| **Location** | `simulation/utils.py` |
| **Parameters** | Per-droplet position lists, time step, droplet count, output path |
| **Returns** | None (writes CSV) |
| **Notes** | Uses forward difference for vx/vy in the CSV (same convention as tracking). |
