# Trajectory Loading

**Module:** `analysis/trajectories.py`

---

## Purpose

Convert a tracking CSV into per-droplet NumPy arrays suitable for ensemble analysis. Apply a full-length filter so every trajectory in the returned dict has the same number of frames.

---

## CSV format expected

```
droplet_id, x, y, vx, vy, theta
0, 320.5, 240.1, 1.2, -0.4, -0.32
0, 321.7, 239.8, ...
1, 150.0, 300.2, ...
```

Column indices `0`, `1`, `2` are read. Columns `3`, `4`, `5` (`vx`, `vy`, `theta`) are ignored — all velocity quantities are recomputed from positions by `analysis.vacf.compute_velocity`.

---

## Full-length filter

Only droplets whose track length equals the global maximum are kept.

```python
max_len = max(v.shape[0] for v in positions_by_id.values())
consistent_ids = [k for k, v in positions_by_id.items() if v.shape[0] == max_len]
```

**Why:** Ensemble averaging over unequal-length arrays requires truncation to the minimum length, which biases the estimate toward the short tracks. The full-length filter gives a clean, unbiased ensemble at the cost of discarding partial tracks.

**Trade-off:** If many droplets appear or disappear mid-video, this filter can retain very few tracks. Use the [diagnostics](diagnostics.md) to assess how much data is discarded.

---

## `load_positions`

```python
def load_positions(csv_path: Path) -> Tuple[PositionsDict, int, int]:
```

**Parameters:**
- `csv_path` — path to a CSV file

**Returns:**
- `positions_consistent` — `Dict[int, ndarray(N, 2)]`, full-length tracks only
- `total_paths` — total number of unique droplet IDs
- `consistent_paths` — number of retained IDs

**Algorithm:**
1. Read all rows, group by `droplet_id`
2. Convert each list of `[x, y]` pairs to `ndarray`
3. Compute `max_len`
4. Return only IDs with `shape[0] == max_len`

---

## `load_positions_with_lengths`

```python
def load_positions_with_lengths(csv_path: Path) -> Tuple[
    PositionsDict, int, int, PositionsDict, Dict[int, int]
]:
```

**Parameters:**
- `csv_path` — path to a CSV file

**Returns:**
- `positions_consistent` — full-length tracks only
- `total_paths` — total unique IDs
- `consistent_paths` — retained IDs count
- `positions_all` — all tracks including partial
- `lengths` — `Dict[int, int]`, frames per droplet ID

**Used by:** `analysis/main.py`, which passes `positions_consistent` to ensemble functions and `positions_all` + `lengths` to diagnostics.

---

## Example usage

```python
from pathlib import Path
from analysis.trajectories import load_positions_with_lengths

csv = Path("output_droplet_tracking/csv_files/video1_droplets.csv")

positions, total, consistent, positions_all, lengths = load_positions_with_lengths(csv)

print(f"Total droplets: {total}")
print(f"Full-length tracks: {consistent}")
print(f"Frames per track: {lengths}")

# positions[0] has shape (N, 2)
r = positions[0]
print(r.shape)  # e.g. (450, 2)
```

---

## Notes

- The function reads the CSV twice (once inside `load_positions` for consistent tracks, once again for `positions_all`). This is acceptable for the file sizes in this project.
- Thread-safety: each call opens and closes the file independently. No global state.
