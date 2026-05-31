# Mean Squared Displacement (MSD)

**Module:** `analysis/msd.py`

---

## Mathematical definition

$$\text{MSD}(\tau) = \left\langle \|\mathbf{r}(t + \tau) - \mathbf{r}(t)\|^2 \right\rangle_t$$

For a 2-D trajectory this expands to:

$$\text{MSD}(\tau) = \left\langle (x(t+\tau)-x(t))^2 + (y(t+\tau)-y(t))^2 \right\rangle_t$$

where `⟨·⟩_t` is a time average over all starting times `t` for which `t + τ ≤ T`.

---

## Theory

### Regime classification

| Behaviour | Exponent α | Physical interpretation |
|-----------|-----------|------------------------|
| Sub-diffusive | α < 1 | Confined, caged, or viscoelastic medium |
| Diffusive | α = 1 | Brownian motion, random walk |
| Super-diffusive | 1 < α < 2 | Active motion with persistence |
| Ballistic | α = 2 | Directed motion (no randomness) |

Active droplets typically show **super-diffusive** behaviour at short to intermediate lag times, crossing over to diffusive at long times after the persistence time τ_p.

### Power-law fit

The dominant analysis is fitting:

$$\text{MSD}(\tau) \sim \tau^\alpha$$

on a log-log scale:

$$\log_{10} \text{MSD} = \alpha \cdot \log_{10} \tau + C$$

α is the slope of this linear fit.

---

## Implementation

### `compute_msd`

```python
def compute_msd(
    r: NDArray[np.floating],
    lag_fraction: float = 1.0,
    minimum_pairs: Optional[int] = None,
) -> NDArray[np.floating]:
```

**Algorithm:**

```python
n = len(r)
max_lag = max(1, int(lag_fraction * (n - 1)))
if minimum_pairs is not None:
    max_lag = min(max_lag, max(1, n - minimum_pairs))

msd = np.zeros(max_lag)
for tau in range(1, max_lag + 1):
    disp = r[tau:] - r[:-tau]          # shape (n-tau, 2)
    msd[tau - 1] = np.mean(np.sum(disp**2, axis=1))
```

- **Estimator:** Direct sliding-window. For lag `τ`, uses all `n - τ` displacement pairs.
- **lag_fraction:** Controls the maximum lag as a fraction of `n-1`. Default `1.0` uses the full range; `0.3` is typical to avoid the noisy long-lag region.
- **minimum_pairs:** If set to e.g. `10`, lags where fewer than 10 pairs exist are excluded. This prevents noisy estimates at long lags.
- **Complexity:** O(N × max_lag). For `lag_fraction=0.3` and `N=500`, this is ≈150,000 operations.

**Returns:** Array of length `max_lag`. Index `i` corresponds to lag `τ = i + 1` (in frame units).

---

### `fit_msd_power_law`

```python
def fit_msd_power_law(
    msd: NDArray[np.floating],
    fit_fraction: float = 0.3,
    minimum_pairs: Optional[int] = None,
) -> Tuple[float, float, float, NDArray[np.floating]]:
```

**Algorithm:**

```python
m = len(msd)
tau = np.arange(1, m + 1)

max_fit = max(1, int(fit_fraction * m))
# limit to [1, max_fit] lag values

log_tau = np.log10(tau_fit[mask])
log_msd = np.log10(msd_fit[mask])

slope, slope_err, r2, fitted_y = _loglog_ols(log_tau, log_msd)
```

- **fit_fraction:** The fraction of the MSD array used for fitting. E.g. `0.3` fits only the first 30% of lag values. This avoids the long-lag regime where the estimator is noisy.
- **OLS:** Uses `statsmodels.OLS` when available (provides standard errors). Falls back to pure-NumPy OLS with identical formulas.
- The fit is on `log10` values, so the slope is the power-law exponent α.

**Returns:** `(alpha, alpha_err, r2, fit_line_full)`

- `fit_line_full` has the same length as `msd`. Values are NaN outside the fit region. This is used by the plotting code to overlay the fit on the full MSD curve.

---

### `ensemble_msd_analysis`

```python
def ensemble_msd_analysis(
    positions_dict: PositionsDict,
    lag_fraction: float = 0.3,
    fit_fraction: float = 0.3,
    minimum_pairs: Optional[int] = None,
) -> Tuple[NDArray, NDArray, float, float, float, NDArray]:
```

**Algorithm:**

```python
# 1. Compute MSD for each droplet
all_msds = [compute_msd(r, lag_fraction, minimum_pairs) for r in positions_dict.values()]

# 2. Truncate to minimum length (handles edge cases)
min_len = min(len(m) for m in all_msds)
ensemble_msd = np.mean([m[:min_len] for m in all_msds], axis=0)

# 3. Fit the ensemble MSD
tau = np.arange(1, min_len + 1)
alpha, alpha_err, r2, fit_line_full = fit_msd_power_law(
    ensemble_msd, fit_fraction=fit_fraction, minimum_pairs=minimum_pairs
)
```

**Key design decision:** Step 3 calls `fit_msd_power_law` internally. This is the same function that notebooks and plotting code call directly. The exponent α is therefore guaranteed to be **bit-identical** regardless of which call path is used.

**Returns:** `(tau, ensemble_msd, alpha, alpha_err, r2, fit_line_full)`

- `tau` is in frame units (integers starting at 1)
- Convert to seconds with `tau_seconds = tau * dt`

---

### `compute_local_msd_exponent`

```python
def compute_local_msd_exponent(
    tau: NDArray[np.floating],
    msd: NDArray[np.floating],
) -> NDArray[np.floating]:
```

**Algorithm:**

```python
return np.gradient(np.log10(msd), np.log10(tau))
```

Uses NumPy's central-difference gradient. The local exponent α(τ) at lag τ is:

$$\alpha(\tau) = \frac{d \log_{10} \text{MSD}}{d \log_{10} \tau}$$

**Inputs:** `tau` and `msd` can be in any consistent positive units (frame integers, seconds — it doesn't matter because both are log-transformed and the ratio is unit-agnostic).

**Validation:** The function raises `ValueError` if inputs are non-finite, non-positive, or non-monotone.

**Interpretation:**
- α(τ) → 2 at very short times: ballistic (free propulsion)
- α(τ) → 1 at long times: diffusive crossover
- The crossover lag is related to the persistence time τ_p

---

### `compute_pair_weighted_msd`

```python
def compute_pair_weighted_msd(
    positions_dict: PositionsDict,
    lag_fraction: float = 0.3,
    minimum_pairs: Optional[int] = None,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]:
```

Unlike `ensemble_msd_analysis` which averages per-droplet MSDs (each droplet contributes equally), this function pools all displacement pairs across all droplets and computes the grand mean:

$$\text{MSD}_\text{weighted}(\tau) = \frac{\sum_{i,t} \|\mathbf{r}_i(t+\tau) - \mathbf{r}_i(t)\|^2}{\sum_{i,t} 1}$$

Droplets with more frames contribute more pairs and therefore more weight. See [Diagnostics](diagnostics.md) for how this is used to assess selection bias.

---

## Example usage

```python
from analysis.msd import ensemble_msd_analysis, compute_local_msd_exponent
from analysis.trajectories import load_positions
from pathlib import Path

positions, _, _ = load_positions(Path("data.csv"))

dt = 0.2  # seconds per frame

tau, ens_msd, alpha, alpha_err, r2, fit = ensemble_msd_analysis(
    positions,
    lag_fraction=0.3,
    fit_fraction=0.3,
)

tau_seconds = tau * dt

print(f"α = {alpha:.3f} ± {alpha_err:.3f},  R² = {r2:.4f}")

# Local exponent
d_alpha = compute_local_msd_exponent(tau_seconds, ens_msd)
```

---

## Example output (`_results.txt`)

```
MSD exponent alpha = 1.473 ± 0.021  R2 = 0.9981
```

---

## Interpretation guide

| α value | Likely cause in active droplets |
|---------|-------------------------------|
| α ≈ 2 | Short-time ballistic; lag < τ_p |
| 1 < α < 2 | Persistent active motion |
| α ≈ 1 | Long-time diffusion; lap >> τ_p |
| α < 1 | Confinement or run-and-tumble with long tumble phases |

---

## Common mistakes

### Mistake 1: Using different `fit_fraction` in batch vs plot

**Before the refactor:** `ensemble_msd_analysis` inlined its fit with one set of parameters while `post_analysis/analysis_utils.fit_msd_power_law` used different defaults.

**After the refactor:** Both are the same function. `ensemble_msd_analysis` calls `fit_msd_power_law` internally with the same `fit_fraction` parameter. They cannot diverge.

### Mistake 2: Fitting the full MSD range

MSD estimates become noisy at long lags because fewer displacement pairs are available. Fitting the full range gives a biased α. Use `fit_fraction=0.3` (fit only the first 30%) as the default.

### Mistake 3: Comparing α across datasets with different `lag_fraction`

If one analysis used `lag_fraction=0.3` and another used `lag_fraction=0.5`, the `ensemble_msd` arrays have different lengths and the fit regions differ. Always use the same `lag_fraction` and `fit_fraction` when comparing.

### Mistake 4: Not converting to seconds

`tau` from `ensemble_msd_analysis` is in frame units. Axes labels and time-scaled quantities require `tau_seconds = tau * dt`. The MSD itself is in position² (pixels² or dimensionless²), unchanged by `dt`.
