# Persistence-Time Analysis

**Modules:** `analysis/vacf.py` (velocity persistence), `analysis/orientation.py` (rotational persistence)

---

## Overview

Two persistence times are extracted from the data:

| Symbol | Name | Source | Model |
|--------|------|---------|-------|
| τ_p | Velocity persistence time | VACF | `A exp(-t/τ_p)` |
| τ_r | Rotational persistence time | Orientation correlation | `exp(-t/τ_r)` |

---

## Translational persistence (τ_p)

### Physical meaning

τ_p characterises how long a droplet maintains its speed and direction before randomising. Formally, for an active Brownian particle:

$$C_v(\tau) = \langle v^2 \rangle \exp\!\left(-\frac{\tau}{\tau_p}\right)$$

A large τ_p means the droplet propels itself persistently in a roughly constant direction for a long time before diffusing. The diffusion coefficient in the long-time limit is:

$$D_\text{eff} = \langle v^2 \rangle \cdot \tau_p$$

### Implementation: `fit_vacf_exponential`

```python
def fit_vacf_exponential(
    tau: NDArray[np.floating],
    vacf: NDArray[np.floating],
    max_fit_fraction: Optional[float] = None,
) -> Tuple[float, float, NDArray, NDArray]:
```

**Model:**

```python
def _model(t, A, tau_p):
    return A * np.exp(-t / tau_p)
```

**Algorithm:**

```python
# 1. Find first zero crossing
negative = np.where(vacf <= 0)[0]
stop = negative[0] if len(negative) > 0 else len(vacf)

# 2. Optionally cap the fit range
if max_fit_fraction is not None:
    stop = min(stop, max(1, int(max_fit_fraction * len(vacf))))

tau_fit = tau[:stop]
vacf_fit = vacf[:stop]

# 3. Fit with scipy.optimize.curve_fit
popt, pcov = curve_fit(_model, tau_fit, vacf_fit, p0=[vacf_fit[0], tau_fit.max()/10])
tau_p = popt[1]
tau_p_err = sqrt(pcov[1, 1])
```

**Initial guess:** `p0 = [vacf_fit[0], tau_fit.max()/10]`. The amplitude is initialised to the VACF at τ=0; the timescale is initialised to 1/10 of the fit range.

**Stopping conditions:**
1. First zero crossing — the exponential model is not valid after the VACF first goes negative
2. `max_fit_fraction` — truncate to a fixed fraction of the full VACF length

**Returns:** `(tau_p, tau_p_err, tau_fit, fit_curve)`

Returns NaN with empty arrays when fewer than 5 points are in the fit region.

**Units:** Same as the input `tau`. In `analysis/main.py`, `tau_vacf_seconds = tau_vacf * dt` is passed, so the returned `tau_p` is in seconds.

---

## Rotational persistence (τ_r)

### Physical meaning

τ_r characterises how quickly the droplet's heading direction randomises. For an active Brownian particle with rotational diffusion coefficient D_r:

$$C_\theta(\tau) = \exp\!\left(-\frac{\tau}{\tau_r}\right), \quad \tau_r = \frac{1}{2 D_r}$$

A large τ_r means the droplet's heading stays nearly constant for a long time.

### Implementation: `fit_orientation_persistence`

```python
def fit_orientation_persistence(
    tau: NDArray[np.floating],
    corr: NDArray[np.floating],
    max_fit_fraction: Optional[float] = None,
) -> Tuple[float, float, NDArray, NDArray]:
```

**Model:**

```python
def _model(t, tau_r):
    return np.exp(-t / tau_r)
```

Note the difference from the VACF fit: the orientation model has **no amplitude prefactor** because `C_θ(0) = 1` always.

**Algorithm:** Identical to `fit_vacf_exponential` except:
- Only one parameter (`tau_r`) instead of two (`A`, `tau_p`)
- `p0 = [tau_fit.max() / 10]`

**Returns:** `(tau_r, tau_r_err, tau_fit, fit_curve)`

**Units:** Same as input `tau` (seconds in `analysis/main.py`).

---

## How persistence times are called in `analysis/main.py`

```python
# VACF persistence
tau_vacf, vacf_raw, vacf_norm = ensemble_vacf_analysis(data, dt)
tau_vacf_seconds = tau_vacf * dt          # convert to seconds

tau_p, tau_p_err, tau_p_fit, vacf_fit_curve = fit_vacf_exponential(
    tau_vacf_seconds,                      # seconds
    vacf_norm,                             # normalised VACF
    max_fit_fraction=vacf_max_fit_fraction
)

# Orientation persistence
tau_orient, orient_corr = ensemble_orientation_corr_analysis(data, dt)
tau_orient_seconds = tau_orient * dt      # convert to seconds

tau_r, tau_r_err, tau_r_fit, orient_fit_curve = fit_orientation_persistence(
    tau_orient_seconds,                    # seconds
    orient_corr,
    max_fit_fraction=orient_max_fit_fraction
)
```

Both are written to `<label>_results.txt`:

```
Velocity persistence tau_p  = 0.423 ± 0.018 s
Rotational persistence tau_r = 0.511 ± 0.024 s
```

---

## Assumptions

1. **Exponential decay:** The VACF and orientation correlation are assumed to decay as single exponentials. Multi-exponential decays (multiple timescales) are not handled.

2. **Stationarity:** The model assumes that τ_p and τ_r are constant in time. If the droplet's dynamics change (e.g. due to chemical depletion in the simulation), the fitted values are effective time-averaged estimates.

3. **Sufficient data:** At least 5 points must be in the fit region. With very short trajectories or fast-decaying correlations, `tau_p` or `tau_r` will be NaN.

4. **No confinement:** Confinement produces a non-zero long-lag VACF floor. The fit to `A exp(-t/τ_p)` will underestimate τ_p in confined systems.

---

## Limitations

- `curve_fit` uses the Levenberg-Marquardt algorithm with no bounds by default. If τ_p or τ_r is initialised far from the true value, convergence may fail, returning NaN.
- Standard errors from `pcov` are valid only under the Gaussian noise assumption and may be unreliable when the residuals are non-Gaussian.
- The VACF is estimated from `N-1` velocity samples. Its statistical uncertainty is not propagated into `tau_p_err`.

---

## Common mistakes

### Mistake 1: Passing frame-index tau instead of seconds

`tau_p` inherits the units of the input `tau`. If you pass frame-index lags (integers), `tau_p` will be in frames. Always pass `tau_seconds = tau * dt`.

### Mistake 2: Fitting the normalised VACF but expecting the raw amplitude

`fit_vacf_exponential` accepts either the raw or normalised VACF. In `analysis/main.py`, the normalised `vacf_norm` is passed, so the fit amplitude `A ≈ 1`. If you pass the raw VACF, `A ≈ C_v(0)`.

### Mistake 3: Setting `max_fit_fraction` too large

If `max_fit_fraction=1.0`, the fit includes the noisy long-lag tail and the oscillatory region after the first zero crossing. The first-zero-crossing stop is always applied first, but `max_fit_fraction` should be set conservatively (e.g. 0.25) for noisy data.
