# Orientation Correlation

**Module:** `analysis/orientation.py`

---

## Mathematical definition

$$C_\theta(\tau) = \left\langle \cos\!\left(\theta(t+\tau) - \theta(t)\right) \right\rangle_t$$

where `θ(t)` is the heading angle of the droplet at time `t`, defined as:

$$\theta(t) = \arctan2\!\left(v_y(t),\, v_x(t)\right)$$

after computing velocity by forward difference and unwrapping phase.

`C_θ(0) = ⟨cos(0)⟩ = 1` always.

---

## Theory

For a persistent random walk with rotational diffusion coefficient `D_r`:

$$C_\theta(\tau) = \exp\!\left(-\frac{\tau}{\tau_r}\right)$$

where the rotational persistence time is `τ_r = 1 / (2 D_r)`.

This exponential decay characterises how quickly the droplet "forgets" its direction of motion. A large `τ_r` means the droplet moves in a nearly straight line for long periods.

**Relationship to VACF persistence:** For an active Brownian particle, `τ_r ≈ τ_p`. In practice they can differ if the speed fluctuates independently of the direction.

---

## Implementation

### `compute_heading_angle`

```python
def compute_heading_angle(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]:
```

```python
v = compute_velocity(r, dt)          # (N-1, 2), from analysis.vacf
theta = np.arctan2(v[:, 1], v[:, 0]) # wrapped to (-π, π)
return np.unwrap(theta)              # remove 2π jumps
```

**Why unwrap?** `arctan2` returns values in `(-π, π)`. A droplet turning continuously would show artificial ±2π jumps without unwrapping. `np.unwrap` adds multiples of 2π to eliminate jumps larger than π.

**Output:** Array of shape `(N-1,)` in radians.

**Note:** `compute_velocity` is imported from `analysis.vacf`, so the velocity definition is shared.

---

### `compute_orientation_corr`

```python
def compute_orientation_corr(
    theta: NDArray[np.floating],
) -> NDArray[np.floating]:
```

**Algorithm:**

```python
n = len(theta)
corr = np.zeros(n)
for tau in range(n):
    corr[tau] = np.mean(np.cos(theta[tau:] - theta[:n - tau]))
```

At each lag τ, compute the cosine of the angle difference between all pairs `(θ(t), θ(t+τ))` and average over all starting times `t`.

**Input:** Heading angles `θ(t)` in radians (unwrapped recommended but not required).  
**Output:** Array of length `n`. Index `i` → lag `τ = i`.  
**Complexity:** O(M²) where M = len(theta).

---

### `ensemble_orientation_corr_analysis`

```python
def ensemble_orientation_corr_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]:
```

**Algorithm:**

```python
all_corr = []
for r in positions_dict.values():
    theta = compute_heading_angle(r, dt)
    all_corr.append(compute_orientation_corr(theta))

min_len = min(len(c) for c in all_corr)
ensemble_corr = np.mean([c[:min_len] for c in all_corr], axis=0)
tau = np.arange(min_len)
```

**Returns:** `(tau, ensemble_corr)`

- `tau` — frame-index lags: 0, 1, …, L-1
- `ensemble_corr` — `C_θ(τ)`, equals 1 at τ=0

**Ensemble strategy:** Mean of per-droplet correlations. Each droplet contributes equally.

---

## Example usage

```python
from analysis.orientation import ensemble_orientation_corr_analysis, fit_orientation_persistence
from analysis.trajectories import load_positions
from pathlib import Path

positions, _, _ = load_positions(Path("data.csv"))
dt = 1.0 / 30.0

tau, orient_corr = ensemble_orientation_corr_analysis(positions, dt)
tau_seconds = tau * dt

tau_r, tau_r_err, tau_r_fit, fit_curve = fit_orientation_persistence(
    tau_seconds, orient_corr, max_fit_fraction=0.3
)
print(f"Rotational persistence: τ_r = {tau_r:.3f} ± {tau_r_err:.3f} s")
```

---

## Interpretation

| Feature | Meaning |
|---------|---------|
| `C_θ(0) = 1` | Perfect self-correlation (always) |
| Fast exponential decay | Low persistence; strong angular noise |
| Slow decay | High persistence; nearly straight trajectories |
| Non-exponential decay | Multiple timescales; complex dynamics |
| Oscillations | Circular or helical motion |

---

## Common mistakes

### Mistake 1: Not unwrapping θ

Without unwrapping, a droplet making a full rotation would produce `cos(±2π) = 1` instead of `cos(very_large_angle) ≈ 0`. The correlation would appear artificially high. `compute_heading_angle` always unwraps.

### Mistake 2: Using raw `arctan2` instead of heading from velocity

`θ = arctan2(vy, vx)` based on the velocity (forward difference) gives the instantaneous direction of motion. Using position angles directly would give the angle to the origin, not the heading direction.

### Mistake 3: Fitting past the first zero crossing

The exponential model is only valid while `C_θ > 0`. `fit_orientation_persistence` automatically truncates at the first non-positive value.
