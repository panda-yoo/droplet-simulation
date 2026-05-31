# Tracking Pipeline

**Module:** `tracking/main.py`, `tracking/utils.py`

---

## Overview

The tracking pipeline converts real microscopy or brightfield videos of active droplets into trajectory CSVs that can be fed directly into the analysis pipeline.

```
videos/*.mp4
    ↓  YOLO object detection + tracking
tracking/utils.generate_csvs()
    ↓
output_droplet_tracking/csv_files/<video_stem>_droplets.csv
    ↓
analysis.main.run_post_processing()
    ↓
output_droplet_tracking/plots/
```

---

## Prerequisites

- A trained YOLO model at `weight/best.pt`
- Video files in `videos/` with extensions `.mp4`, `.mov`, `.avi`, or `.mkv`
- Python packages: `ultralytics`

---

## Entry point: `run_tracking`

```python
# tracking/main.py
def run_tracking(
    output_dir: str = "output_droplet_tracking",
    save_annotated: bool = True,
) -> None:
```

**Parameters:**
- `output_dir` — root directory for all outputs
- `save_annotated` — if True, saves YOLO-annotated videos with bounding boxes

**What it does:**

```python
model = YOLO(weight)                  # load weight/best.pt

video_paths = sorted([...])           # discover *.mp4, *.mov, *.avi, *.mkv

generate_csvs(
    model=model,
    video_paths=video_paths,
    output_dir=output_path,
    save_annotated=save_annotated,
)

run_post_processing(
    base_dir=str(csv_dir),            # output_droplet_tracking/csv_files/
    output_dir=str(plots_dir),        # output_droplet_tracking/plots/
    dt=dt,                            # 1.0 / 30.0 (hardcoded FPS = 30)
)
```

**Time step:** `dt = 1.0 / fps` where `fps = 30.0` is hardcoded in `tracking/main.py`. If your videos were recorded at a different frame rate, modify `fps` in that file.

---

## `generate_csvs`

```python
def generate_csvs(
    model: YOLO,
    videos_dir: Path,
    video_exts: Set[str],
    video_paths: List[Path],
    output_dir: Path,
    save_annotated: bool = True,
) -> None:
```

**For each video:**

1. **Run YOLO tracking:**
   ```python
   track_results = model.track(source=str(video_path), stream=True, ...)
   ```
   YOLO's `.track()` applies ByteTrack/BotSort to associate detections across frames.

2. **Collect bounding-box centres:**
   ```python
   for r in track_results:
       ids = r.boxes.id.cpu().numpy()
       xyxy = r.boxes.xyxy.cpu().numpy()
       for obj_id, box in zip(ids, xyxy):
           xc = (box[0] + box[2]) / 2
           yc = (box[1] + box[3]) / 2
           positions_by_id[int(obj_id)].append([xc, yc])
   ```
   Position units are pixels (image coordinates).

3. **Write CSV:**
   ```python
   rows = positions_to_rows(positions_by_id)
   ```
   See below.

4. **Move annotated video** (if `save_annotated=True`):
   YOLO writes to a temporary `temp/` directory; the file is renamed and moved to `annotated_videos/`.

5. **Clean up temporary YOLO directory.**

---

## `positions_to_rows`

```python
def positions_to_rows(
    positions_by_id: Dict[int, List[List[float]]],
    dt: float = 1.0,
) -> List[List[float]]:
```

For each droplet ID:
1. Convert positions to `ndarray(N, 2)`
2. Compute velocity for CSV columns only:
   ```python
   v = _velocity_for_csv(pos, dt)   # uses compute_velocity + last-row padding
   theta = np.arctan2(v[:, 1], v[:, 0])
   ```
3. Emit one row per frame: `[droplet_id, x, y, vx, vy, theta]`

**Important:** The `vx`, `vy`, `theta` columns written here are **not read back** by `analysis/trajectories.py`. Analysis always recomputes velocity from `(x, y)` using `compute_velocity` (forward difference, shape N-1). The stored columns are informational only (e.g. for quick inspection of the CSV).

---

## `_velocity_for_csv`

```python
def _velocity_for_csv(pos, dt=1.0):
    v_inner = compute_velocity(pos, dt)    # shape (N-1, 2)
    return np.vstack([v_inner, v_inner[-1:]])  # shape (N, 2), last row duplicated
```

This pads the velocity array to the same length as positions so the CSV has equal row counts. The last velocity value is duplicated for the final frame because there is no subsequent frame to difference with.

---

## Output directories

```
output_droplet_tracking/
    csv_files/
        video1_droplets.csv
        video2_droplets.csv
    annotated_videos/
        annotated_video1.mp4
        annotated_video2.mp4
    plots/
        video1_loglog_msd.png
        video1_results.txt
        ...
```

---

## Modifying the FPS

`fps` is hardcoded in `tracking/main.py`:

```python
fps: float = 30.0
dt: float = 1.0 / fps
```

To use a different frame rate, change `fps`. Alternatively, use `DropletPipeline.run_analysis()` with a custom `dt` to reanalyse existing CSVs without re-running tracking.

---

## Using `run_tracking` from a script

```python
from tracking.main import run_tracking

run_tracking(
    output_dir="my_output",
    save_annotated=True,
)
```

Or via the pipeline:

```python
from main import DropletPipeline
p = DropletPipeline()
p.run_tracking(output_dir="my_output", save_annotated=True)
```
