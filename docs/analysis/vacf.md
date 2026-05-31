# Velocity Autocorrelation Function (VACF)

**Module:** `analysis/vacf.py`

---

## Mathematical definition

$$C_v(\tau) = \left\langle \mathbf{v}(t) \cdot \mathbf{v}(t + \tau) \right\rangle_t$$

Expanding the dot product for a 2-D trajectory:

$$C_v(\tau) = \left\langle v_x(t) v_x(t+\tau) + v_y(t) v_y(t+\tau) \right\rangle_t$$

where `⟨·⟩_t` is a time average over all starting times `t` for which `t + τ ≤ T_v` (`T_v` = total velocity frames = `N - 1`).

At `τ = 0`:

$$C_v(0) = \left\langle v_x^2 + v_y^2 \right\rangle = \left\langle |\mathbf{v}|^2 \right\rangle$$

This is proportional to the mean kinetic energy per unit mass.

The normalised VACF is:

$$\hat{C}_v(\tau) = \frac{C_v(\tau)}{C_v(0)}$$

which equals 1 at τ=0 and decays toward 0.

---

## Velocity convention

**Canonical definition** (used everywhere in this project):

```python
# analysis/vacf.py :: compute_velocity()
v = (r[1:] - r[:-1]) / dt      # shape (N-1, 2)
```

This is a **forward difference** at each frame interval. The result has shape `(N-1, 2)` — one fewer row than the position array.

**This is the only velocity definition used by the analysis package.** It is imported by `orientation.py`, `psd.py`, `diagnostics.py`, and `tracking/utils.py`.

The `vx/vy/theta` columns stored in the CSV are computed with a padded version of this formula (last row duplicated to maintain equal row count) and are **not read back** during analysis.

---

## Implementation

### `compute_velocity`

```python
def compute_velocity(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]:
```

```python
return (r[1:] - r[:-1]) / dt
```

**Input:** `(N, 2)` position array  
**Output:** `(N-1, 2)` velocity array  
**Complexity:** O(N)

---

### `compute_vacf`

```python
def compute_vacf(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]:
```

**Algorithm:**

```python
v = compute_velocity(r, dt)      # shape (N-1, 2)
n = len(v)
vacf = np.zeros(n)

for tau in range(n):
    corr = (v[:n-tau, 0] * v[tau:, 0]
          + v[:n-tau, 1] * v[tau:, 1])
    vacf[tau] = np.mean(corr)
```

- Uses all `n - τ` velocity pairs for each lag
- Returns array of length `n = N - 1`
- Index `i` corresponds to lag `τ = i` (in frame units, starting at 0)
- `vacf[0] = C_v(0) = ⟨|v|²⟩`

**Complexity:** O((N-1)²). For `N = 500` this is ~250,000 operations.

---

### `ensemble_vacf_analysis`

```python
def ensemble_vacf_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray[np.int_], NDArray[np.floating], NDArray[np.floating]]:
```

**Algorithm:**

```python
all_vacf = [compute_vacf(r, dt) for r in positions_dict.values()]
min_len = min(len(v) for v in all_vacf)
ensemble_vacf = np.mean([v[:min_len] for v in all_vacf], axis=0)
tau = np.arange(min_len)

vacf_norm = ensemble_vacf / ensemble_vacf[0] if ensemble_vacf[0] != 0 else ensemble_vacf
```

**Returns:** `(tau, vacf_raw, vacf_norm)`

- `tau` — frame-index lags: 0, 1, …, L-1
- `vacf_raw` — raw VACF in velocity² units
- `vacf_norm` — divided by `C_v(0)`; equals 1 at τ=0

**Ensemble strategy:** Mean of per-droplet VACFs. Each droplet contributes equally regardless of its speed.

---

## Example usage

```python
from analysis.vacf import ensemble_vacf_analysis, fit_vacf_exponential
from analysis.trajectories import load_positions
from pathlib import Path

positions, _, _ = load_positions(Path("data.csv"))
dt = 1.0 / 30.0  # 30 fps video

tau, vacf_raw, vacf_norm = ensemble_vacf_analysis(positions, dt)
tau_seconds = tau * dt

# Fit persistence time
tau_p, tau_p_err, tau_p_fit, fit_curve = fit_vacf_exponential(
    tau_seconds, vacf_norm, max_fit_fraction=0.25
)
print(f"Velocity persistence time: τ_p = {tau_p:.3f} ± {tau_p_err:.3f} s")
```

---

## Interpretation

| VACF feature | Physical meaning |
|-------------|-----------------|
| Rapid decay to 0 | Short persistence; nearly diffusive |
| Slow exponential decay | Long persistence; strong self-propulsion |
| Oscillation / negative lobe | Velocity reversal; run-and-tumble or confinement |
| `C_v(0)` value | Mean speed squared |

For a simple persistent random walk, the VACF decays as `exp(-t/τ_p)` where τ_p is the persistence time (see [Persistence](persistence.md)).

---

## Relationship to MSD

For a persistent random walk, the Green-Kubo relation links VACF and MSD:

$$\text{MSD}(\tau) = 2 \int_0^\tau (\tau - s) C_v(s) \, ds$$

In practice this is not computed directly, but the relationship means the persistence time extracted from the VACF should be consistent with the crossover in the local MSD exponent α(τ).

---

## Common mistakes

### Mistake 1: Using the wrong velocity definition

If velocity is computed differently in different places (e.g. central difference vs forward difference), the VACF will differ. In this codebase `compute_velocity` is imported everywhere — it cannot differ.

### Mistake 2: Plotting the full VACF length

The VACF at long lags is estimated from very few velocity pairs and is very noisy. `analysis/main.py` plots only the first quarter: `cut = len(tau_vacf) // 4`.

### Mistake 3: Fitting the VACF after the first zero crossing

The exponential model `A exp(-t/τ_p)` is only valid in the positive-valued regime. `fit_vacf_exponential` automatically stops at the first non-positive value.
