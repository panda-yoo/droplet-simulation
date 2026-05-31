# Troubleshooting Guide

Debugging workflows for common problems in the active-droplet analysis pipeline.

---

## MSD exponent mismatch

**Symptom:** The α printed during batch analysis differs from the α shown in a notebook plot.

**Cause before the refactor:** `ensemble_msd_analysis` computed its own inline OLS fit, while `post_analysis/analysis_utils.fit_msd_power_law` was a separate function. Any difference in `fit_fraction` or the fit region would produce different α values.

**Current status (after refactor):** `ensemble_msd_analysis` calls `fit_msd_power_law` internally. Both code paths produce bit-identical α.

**Debugging workflow:**

```python
from analysis.msd import ensemble_msd_analysis, fit_msd_power_law
import numpy as np

# positions loaded from the same CSV in both paths
tau, ens_msd, alpha_batch, *_ = ensemble_msd_analysis(
    positions, lag_fraction=0.3, fit_fraction=0.3
)
alpha_plot, *_ = fit_msd_power_law(ens_msd, fit_fraction=0.3)

print(abs(alpha_batch - alpha_plot))   # must be < 1e-12
```

**If they still differ:**
1. Are you passing the same `fit_fraction` to both calls?
2. Is the notebook using an old `post_analysis/analysis_utils.py` that was not yet updated? Restart the kernel and reimport.
3. Check `sys.modules` for stale imports: `print(sys.modules.get('analysis.msd'))`.

---

## VACF mismatch

**Symptom:** VACF shape or values differ between a notebook and the batch pipeline.

**Debugging workflow:**

```python
from analysis.vacf import compute_velocity, ensemble_vacf_analysis

# Check velocity shape
r = positions[0]
v = compute_velocity(r, dt=dt)
print(f"r.shape={r.shape}, v.shape={v.shape}")  # should be (N,2) and (N-1,2)

# Check VACF
tau, vacf_raw, vacf_norm = ensemble_vacf_analysis(positions, dt=dt)
print(f"vacf[0]={vacf_raw[0]:.4f}  (should be mean squared speed)")
print(f"vacf_norm[0]={vacf_norm[0]:.4f}  (should be 1.0)")
```

**Common causes:**
1. **dt mismatch:** If `dt` differs between calls, `compute_velocity` will scale velocities differently, changing `C_v(0)` and τ_p. Always pass the same `dt`.
2. **Old shim import:** If `from analysis.compute import compute_vacf` is used in a notebook after a kernel restart, it should resolve to the canonical function via the shim. If not, there may be a stale `.pyc` file — delete `__pycache__/` and restart.

---

## PSD mismatch

**Symptom:** β differs between runs, or `fit_psd_powerlaw` returns NaN.

**Debugging workflow:**

```python
from analysis.psd import ensemble_psd_analysis, fit_psd_powerlaw

freq, psd_vx, psd_vy = ensemble_psd_analysis(positions, dt=dt)
print(f"freq range: [{freq[0]:.3f}, {freq[-1]:.3f}] Hz")
print(f"Nyquist = {0.5/dt:.3f} Hz")

beta, beta_err, r2, fit_freqs, fit_line = fit_psd_powerlaw(
    freq, psd_vx, fmin=0.1, fmax=2.0
)
print(f"Points in fit range: {len(fit_freqs)}")
print(f"beta = {beta}")
```

**Common causes:**
1. **No points in [fmin, fmax]:** If `fmax > 0.5/dt` (above Nyquist) or if the frequency resolution is too coarse, `len(fit_freqs) < 3` → NaN returned. Fix by checking `freq` range and adjusting `fmin`/`fmax`.
2. **Wrong dt:** If `dt` is 10× too small, all frequencies are 10× too high and the fit window misses the correct region.
3. **Too-short trajectory:** If `len(v) < 8`, `compute_psd_components` returns empty arrays. Check that your trajectories have at least 10 frames.

---

## Inconsistent `dt`

**Symptom:** All time-scaled quantities (τ_p, τ_r, PSD frequencies, tau_seconds axis) are wrong by a constant factor.

**Diagnosis:**

```python
# Check what dt the pipeline used:
# In analysis/main.py::post_processing, dt is a parameter.
# In tracking/main.py: dt = 1.0 / 30.0
# In simulation/main.py: dt_val = 0.002

# Verify:
fps = 30.0
dt_expected = 1.0 / fps
print(f"Expected dt = {dt_expected:.6f} s per frame")
```

**Fix:** Always pass the correct `dt` matching how the video was recorded or the simulation was run. MSD exponent α is unaffected because `tau_seconds = tau * dt` appears on both log axes and `dt` cancels in the slope.

---

## `lag_fraction` issues

**Symptom:** MSD curves are truncated unexpectedly or the fit region is empty.

**Diagnosis:**

```python
r = positions[0]
n = len(r)
lag_fraction = 0.3
max_lag = max(1, int(lag_fraction * (n - 1)))
print(f"n={n}, max_lag={max_lag}")

from analysis.msd import compute_msd
msd = compute_msd(r, lag_fraction=lag_fraction)
print(f"msd.shape={msd.shape}")  # should be (max_lag,)
```

**Common causes:**
1. `lag_fraction=1.0` with short trajectories (`n < 20`) produces very noisy long-lag estimates.
2. Different `lag_fraction` values in different parts of the code. Both `ensemble_msd_analysis` and `compute_pair_weighted_msd` take `lag_fraction` as a parameter — ensure they use the same value when comparing results.

---

## `fit_fraction` issues

**Symptom:** α is too high or too low; fit appears to include the noisy tail.

**Diagnosis:**

```python
from analysis.msd import fit_msd_power_law
import matplotlib.pyplot as plt

msd = ...   # your ensemble MSD
tau = np.arange(1, len(msd) + 1)

for ff in [0.1, 0.2, 0.3, 0.5, 1.0]:
    alpha, _, _, fit_line = fit_msd_power_law(msd, fit_fraction=ff)
    fit_mask = ~np.isnan(fit_line)
    plt.loglog(tau[fit_mask], msd[fit_mask])
    print(f"fit_fraction={ff}: alpha={alpha:.3f}, fit points={fit_mask.sum()}")

plt.show()
```

**Recommendation:** Use `fit_fraction=0.3` as the default. Increase it only if the trajectory is very long and the MSD is clean at long lags. Decrease it if the MSD bends away from power-law at intermediate lags.

---

## Ensemble averaging differences

**Symptom:** Ensemble MSD with N droplets differs significantly from a single-droplet MSD.

**This is expected** — the ensemble average reduces noise. However, if the ensemble MSD exponent α differs dramatically from individual droplet exponents, investigate:

```python
from analysis.msd import compute_msd, fit_msd_power_law
import numpy as np

individual_alphas = []
for r in positions.values():
    msd = compute_msd(r, lag_fraction=0.3)
    alpha, *_ = fit_msd_power_law(msd, fit_fraction=0.3)
    individual_alphas.append(alpha)

print(f"Individual alpha: mean={np.mean(individual_alphas):.3f}, "
      f"std={np.std(individual_alphas):.3f}")
```

If `std` is large, individual droplets have very different behaviour and the ensemble average may not be representative.

---

## Track filtering problems

**Symptom:** Very few tracks are retained (e.g. only 1–2 out of 20).

**Diagnosis:**

```python
from analysis.trajectories import load_positions_with_lengths
from analysis.diagnostics import retained_vs_discarded_lengths
from pathlib import Path

_, total, consistent, _, lengths = load_positions_with_lengths(Path("data.csv"))
ret_len, disc_len, n_ret, n_disc = retained_vs_discarded_lengths(lengths)

print(f"Total: {total},  Retained: {n_ret},  Discarded: {n_disc}")
print(f"Max length: {max(lengths.values())}")
print(f"Length distribution: {sorted(set(lengths.values()))}")
```

**Common causes:**
1. **One very long track:** If one droplet is tracked for all 500 frames but the others disappear at frame 400, only that one droplet is retained. In this case `n_ret=1`.
2. **YOLO ID switching:** If YOLO reassigns track IDs mid-video, one physical droplet may appear as multiple short tracks, all discarded.
3. **Fix options:**
   - Use `pair_weighted_msd` which includes all tracks
   - Use `positions_all` directly and set a minimum length threshold manually

---

## Stale `__pycache__`

**Symptom:** Changes to source files are not reflected in notebook output.

**Fix:**

```bash
find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null
find . -name "*.pyc" -delete
```

Then restart the Jupyter kernel.

---

## `statsmodels` not installed

**Symptom:** `alpha_err` and `r2` are computed with the NumPy fallback.

The NumPy fallback produces identical `alpha`/`beta` point estimates. The `alpha_err` from the fallback uses the same OLS formula as `statsmodels` and is equivalent under the Gaussian noise assumption.

**To verify which path is used:**

```python
import statsmodels.api as sm   # if this raises ImportError, fallback is used
```

**To install:**

```bash
pip install statsmodels
```

---

## `scipy.optimize.curve_fit` failure (NaN tau_p or tau_r)

**Symptom:** `_results.txt` shows `tau_p = nan`.

**Causes:**
1. VACF goes negative immediately (at τ=1 frame). Only `tau_p_fit` with 0 or 1 points → NaN.
2. Initial guess too far from truth. The initial `tau_p` guess is `tau_fit.max() / 10`.

**Debugging:**

```python
from analysis.vacf import ensemble_vacf_analysis, fit_vacf_exponential
import numpy as np

tau, vacf_raw, vacf_norm = ensemble_vacf_analysis(positions, dt=dt)
tau_seconds = tau * dt

# Print first zero crossing
negative = np.where(vacf_norm <= 0)[0]
print(f"First zero crossing at tau index: {negative[0] if len(negative)>0 else 'none'}")
print(f"Fit will use {negative[0] if len(negative)>0 else len(vacf_norm)} points")
```

If the zero crossing is at index 1 or 2, the VACF has essentially no positive regime. This means there is no persistent motion — τ_p is not well-defined for this data.
