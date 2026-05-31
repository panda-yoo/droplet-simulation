# Power Spectral Density (PSD)

**Module:** `analysis/psd.py`

---

## Mathematical definition

The power spectral density of velocity component `v_x(t)` is:

$$S_{v_x}(f) = \lim_{T \to \infty} \frac{1}{T} \left| \int_0^T v_x(t) e^{-2\pi i f t} \, dt \right|^2$$

In practice it is estimated from the finite discrete-time series using Welch's method.

The power-law model fitted to the PSD is:

$$S(f) \sim f^{-\beta}$$

where β > 0 means the PSD decreases with frequency (reddish noise). β = 0 is white noise. β = 2 corresponds to 1/f² (Brownian noise from a Langevin process).

---

## Theory

For an active Brownian particle with exponential velocity correlations (VACF ~ exp(-t/τ_p)), the velocity PSD is a Lorentzian:

$$S_v(f) = \frac{2 \langle v^2 \rangle \tau_p}{1 + (2\pi f \tau_p)^2}$$

At low frequencies (`f ≪ 1/τ_p`): `S ≈ const` (white)  
At high frequencies (`f ≫ 1/τ_p`): `S ~ f^{-2}`

The fitted β therefore depends on which frequency range is fitted (`fmin`, `fmax`).

---

## Implementation

### `compute_psd_components`

```python
def compute_psd_components(
    r: NDArray[np.floating],
    dt: float = 1.0,
) -> Tuple[NDArray, NDArray, NDArray]:
```

**Algorithm:**

```python
v = compute_velocity(r, dt)          # (N-1, 2), from analysis.vacf
vx = v[:, 0] - np.mean(v[:, 0])     # mean-subtract
vy = v[:, 1] - np.mean(v[:, 1])

fs = 1.0 / dt                        # sampling frequency (Hz)
nperseg = min(256, len(vx))          # Welch segment length

freq, psd_vx = welch(vx, fs=fs, nperseg=nperseg)
_, psd_vy    = welch(vy, fs=fs, nperseg=nperseg)
```

**Why mean-subtract?** The DC component at f=0 would otherwise dominate and distort the PSD. Mean subtraction removes the time-averaged drift.

**Welch's method:** Divides the signal into overlapping segments, computes a periodogram for each, and averages. This reduces the variance of the PSD estimate at the cost of frequency resolution. Segment length `nperseg = min(256, len(vx))` balances resolution and variance.

**Frequency units:** `freq` is in cycles per unit of `dt`. If `dt` is in seconds then `freq` is in Hz. This convention is enforced throughout — no manual rescaling is needed.

**Returns:** `(freq, psd_vx, psd_vy)`  
Returns three empty arrays if `len(vx) < 8` (minimum for Welch to be meaningful).

---

### `fit_psd_powerlaw`

```python
def fit_psd_powerlaw(
    freqs: NDArray[np.floating],
    psd: NDArray[np.floating],
    fmin: float = 0.1,
    fmax: float = 2.0,
) -> Tuple[float, float, float, NDArray, NDArray]:
```

**Algorithm:**

```python
mask = (freqs >= fmin) & (freqs <= fmax) & (psd > 0)
logf = np.log10(freqs[mask])
logp = np.log10(psd[mask])

slope, slope_err, r2, fit_line = _loglog_ols(logf, logp)
beta = -slope    # PSD ~ f^{-beta}, so slope of log-log = -beta
```

**fit_fraction vs fmin/fmax:** Unlike MSD fitting which uses a fraction of the lag range, PSD fitting uses explicit frequency bounds `[fmin, fmax]`. Default `fmin=0.1 Hz`, `fmax=2.0 Hz`. You must set these to the physically meaningful power-law regime for your system.

**Returns:** `(beta, beta_err, r2, fit_freqs, fit_line)`

- `beta` — positive exponent (PSD decays as f^{−β})
- `fit_freqs` — the frequency values within [fmin, fmax] used for fitting
- `fit_line` — fitted PSD values at `fit_freqs`

Returns NaN and empty arrays if fewer than 3 points are in the fit range.

---

### `ensemble_psd_analysis`

```python
def ensemble_psd_analysis(
    positions_dict: PositionsDict,
    dt: float = 1.0,
) -> Tuple[NDArray, NDArray, NDArray]:
```

**Algorithm:**

```python
all_psd_vx, all_psd_vy = [], []
for r in positions_dict.values():
    freq, psd_vx, psd_vy = compute_psd_components(r, dt)
    if len(freq) == 0: continue
    all_psd_vx.append(psd_vx)
    all_psd_vy.append(psd_vy)

min_len = min(len(p) for p in all_psd_vx)
ensemble_psd_vx = np.mean([p[:min_len] for p in all_psd_vx], axis=0)
ensemble_psd_vy = np.mean([p[:min_len] for p in all_psd_vy], axis=0)
```

Trajectories with fewer than 8 velocity samples are skipped. All PSD arrays are truncated to the shortest frequency axis before averaging.

**Returns:** `(freq, ensemble_psd_vx, ensemble_psd_vy)`

---

## Example usage

```python
from analysis.psd import ensemble_psd_analysis, fit_psd_powerlaw
from analysis.trajectories import load_positions
from pathlib import Path

positions, _, _ = load_positions(Path("data.csv"))
dt = 1.0 / 30.0

freq, psd_vx, psd_vy = ensemble_psd_analysis(positions, dt)

beta, beta_err, r2, fit_freqs, fit_line = fit_psd_powerlaw(
    freq, psd_vx, fmin=0.1, fmax=5.0
)
print(f"PSD exponent β = {beta:.3f} ± {beta_err:.3f}")
```

---

## Interpretation

| β value | Physical meaning |
|---------|-----------------|
| β ≈ 0 | White noise velocity fluctuations |
| β ≈ 1 | 1/f noise (long-memory correlations) |
| β ≈ 2 | Brownian (Ornstein-Uhlenbeck) velocity process |
| β > 2 | Very strong low-frequency dominance |

---

## Common mistakes

### Mistake 1: Wrong dt (wrong frequency axis)

`freq` in the output is `[0, 1/(N·dt), 2/(N·dt), ...]` in Hz. If `dt` is wrong, all frequency values are wrong. `fmin` and `fmax` in `fit_psd_powerlaw` are in the same units as `freq`, so they must match the actual `dt` used.

### Mistake 2: Fitting across the Nyquist frequency

The maximum valid frequency is `f_Nyquist = 1/(2·dt)`. Setting `fmax` above this will include aliased content. Check that `fmax < 1/(2·dt)`.

### Mistake 3: Comparing β across datasets with different segment lengths

If one dataset has `len(vx)=100` and another has `len(vx)=1000`, they use different `nperseg` values (100 vs 256), giving different frequency resolutions. The fit range [fmin, fmax] may span different numbers of points, affecting the fit.
