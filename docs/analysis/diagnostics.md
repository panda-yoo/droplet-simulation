# Ensemble Bias Diagnostics

**Module:** `analysis/diagnostics.py`

---

## Motivation

The analysis pipeline retains only **full-length tracks** — trajectories whose length equals the global maximum in a CSV file. This filter is necessary for clean ensemble averaging but introduces potential selection bias: if longer-lived droplets are faster, larger, or more persistent than shorter-lived ones, the retained ensemble is not representative of all droplets.

The four diagnostic functions quantify this bias **without changing the analysis results**.

---

## `retained_vs_discarded_lengths`

```python
def retained_vs_discarded_lengths(
    lengths: Dict[int, int],
) -> Tuple[List[int], List[int], int, int]:
```

**Parameters:**
- `lengths` — mapping `droplet_id → track_length_in_frames` (from `load_positions_with_lengths`)

**Returns:** `(retained_lengths, discarded_lengths, n_retained, n_discarded)`

**Algorithm:**

```python
max_len = max(lengths.values())
retained  = [l for l in lengths.values() if l == max_len]
discarded = [l for l in lengths.values() if l != max_len]
```

A track is "retained" if its length equals `max_len`. A track is "discarded" if it is shorter.

**Output:** Used in `analysis/main.py` to:
1. Print retained/discarded counts to stdout
2. Pass to `plot_track_length_hist` (histogram of retained vs discarded lengths)
3. Write counts to `<label>_results.txt`

**Interpretation:** If `n_discarded >> n_retained`, the full-length filter is very aggressive. Consider whether partial tracks might be worth including with a different strategy (e.g. pair-weighted MSD).

---

## `retained_vs_discarded_speed_stats`

```python
def retained_vs_discarded_speed_stats(
    positions_all: PositionsDict,
    lengths: Dict[int, int],
    dt: float = 1.0,
) -> Dict[str, float]:
```

**Parameters:**
- `positions_all` — all trajectories including partial tracks (from `load_positions_with_lengths`)
- `lengths` — track lengths per droplet ID
- `dt` — time step between frames

**Returns:** dict with keys:
- `retained_mean_speed` — mean of per-droplet time-averaged speeds (retained tracks)
- `retained_std_speed` — std deviation of those speeds
- `discarded_mean_speed` — same for discarded tracks
- `discarded_std_speed`
- `n_retained`
- `n_discarded`

**Algorithm:**

```python
for did, r in positions_all.items():
    v = compute_velocity(r, dt)
    mean_speed = np.mean(np.sqrt(v[:,0]**2 + v[:,1]**2))
    if lengths[did] == max_len:
        retained_speeds.append(mean_speed)
    else:
        discarded_speeds.append(mean_speed)
```

**Interpretation:**
- If `retained_mean_speed ≈ discarded_mean_speed` → no speed-based selection bias.
- If `retained_mean_speed >> discarded_mean_speed` → the retained ensemble contains only fast droplets, which will overestimate persistence and diffusion coefficients.
- Written to `<label>_results.txt` under `--- Bias diagnostics ---`.

---

## `pair_weighted_msd`

This is a re-export of `analysis.msd.compute_pair_weighted_msd`:

```python
pair_weighted_msd = compute_pair_weighted_msd
```

See [MSD documentation](msd.md#compute_pair_weighted_msd) for the implementation.

**Purpose as diagnostic:** Compare the pair-weighted MSD (all tracks, weighted by number of displacement pairs) against the standard ensemble MSD (full-length tracks only, each droplet weighted equally). A large difference suggests selection bias.

**In `analysis/main.py`:**

```python
tau_weighted, msd_weighted = compute_pair_weighted_msd(
    positions_all,           # includes partial tracks
    lag_fraction=msd_lag_fraction,
    minimum_pairs=msd_minimum_pairs,
)

# Plotted alongside the standard ensemble MSD
plot_msd_comparison(
    tau_seconds, ensemble_msd,            # from full-length tracks only
    tau_weighted_seconds, msd_weighted,   # from all tracks
    ...
)
```

**Saved as:** `<label>_msd_pair_weighted.png`

---

## `pair_weighted_vacf`

```python
def pair_weighted_vacf(
    positions_dict: PositionsDict,
    dt: float = 1.0,
    lag_fraction: float = 0.3,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]:
```

**Parameters:**
- `positions_dict` — any collection of trajectories (pass `positions_all` for bias diagnosis)
- `dt` — time step
- `lag_fraction` — fraction of shortest track used as maximum lag

**Algorithm:**

```python
for r in positions_dict.values():
    v = compute_velocity(r, dt)
    for tau in range(max_lag + 1):
        dot = v[:nv-tau, 0] * v[tau:, 0] + v[:nv-tau, 1] * v[tau:, 1]
        corr_sum[tau] += np.sum(dot)
        pair_counts[tau] += len(dot)

vacf_weighted = corr_sum / pair_counts
```

**Returns:** `(tau, vacf_weighted)` — lag in frame units, pair-weighted VACF

**Purpose:** Like `pair_weighted_msd`, compare against `ensemble_vacf_analysis` to detect whether the full-length track filter biases the velocity correlation estimate.

**Note:** `pair_weighted_vacf` is not called by `analysis/main.py` by default. Import and call it directly in a notebook when investigating bias.

---

## Complete diagnostic workflow example

```python
from analysis.trajectories import load_positions_with_lengths
from analysis.msd import ensemble_msd_analysis
from analysis.diagnostics import (
    retained_vs_discarded_lengths,
    retained_vs_discarded_speed_stats,
    pair_weighted_msd,
    pair_weighted_vacf,
)
from pathlib import Path

positions, total, consistent, positions_all, lengths = load_positions_with_lengths(
    Path("data.csv")
)

dt = 0.2

# 1. How many tracks were discarded?
ret_len, disc_len, n_ret, n_disc = retained_vs_discarded_lengths(lengths)
print(f"Retained: {n_ret}  Discarded: {n_disc}")

# 2. Are retained tracks faster than discarded ones?
stats = retained_vs_discarded_speed_stats(positions_all, lengths, dt)
print(f"Retained mean speed : {stats['retained_mean_speed']:.4f}")
print(f"Discarded mean speed: {stats['discarded_mean_speed']:.4f}")

# 3. Does pair-weighted MSD differ from ensemble MSD?
tau, ens_msd, alpha, *_ = ensemble_msd_analysis(positions, lag_fraction=0.3, fit_fraction=0.3)
tau_w, msd_w = pair_weighted_msd(positions_all, lag_fraction=0.3)

import matplotlib.pyplot as plt
plt.figure()
plt.loglog(tau * dt, ens_msd, label=f"Ensemble (α={alpha:.2f})")
plt.loglog(tau_w * dt, msd_w, label="Pair-weighted (all tracks)")
plt.legend(); plt.show()
```

---

## Interpreting the bias

| Diagnostic | Low bias | High bias |
|------------|----------|-----------|
| `n_discarded / n_retained` | < 0.5 | > 2 |
| `retained_mean_speed - discarded_mean_speed` | ≈ 0 | > 20% of retained mean |
| pair-weighted MSD vs ensemble MSD | overlap within noise | systematically different slope |

If high bias is detected, consider:
1. Using pair-weighted MSD/VACF instead of ensemble averages for reporting
2. Relaxing the full-length filter (e.g. keep tracks ≥ 80% of max length)
3. Reporting the bias as a systematic uncertainty in the paper
