"""
tracking.utils
==============
Helper functions for YOLO-based droplet tracking.

Contents
--------
- positions_to_rows : convert per-droplet positions → flat CSV rows
- generate_csvs     : run YOLO tracking on videos and write CSVs

Notes
-----
Velocity in the CSV is a padded forward-difference used only for the
vx/vy/theta columns stored in the file.  All analysis quantities (MSD,
VACF, PSD, orientation correlation) are recomputed from the (x, y) position
columns by ``analysis.vacf.compute_velocity`` — the sole canonical velocity
definition used throughout the project.
"""

import csv
import numpy as np
from numpy.typing import NDArray
from pathlib import Path
from typing import Dict, List, Set

from ultralytics import YOLO

from analysis.vacf import compute_velocity


# =========================================================
# CSV velocity helper (padded to match position length)
# =========================================================

def _velocity_for_csv(
    pos: NDArray[np.floating],
    dt: float = 1.0,
) -> NDArray[np.floating]:
    """
    Return forward-difference velocity padded to the same length as ``pos``.

    Parameters
    ----------
    pos : ndarray, shape (N, 2)
        Ordered (x, y) positions.
    dt : float
        Time step between consecutive frames.

    Returns
    -------
    v : ndarray, shape (N, 2)
        Forward-difference velocity; the last row is a copy of the second-to-last
        so that the CSV has the same number of rows as the position array.

    Notes
    -----
    This function is used ONLY for writing the vx/vy/theta columns to the CSV.
    Analysis code uses ``analysis.vacf.compute_velocity`` directly on positions.
    """
    if pos.shape[0] < 2:
        return np.zeros_like(pos)
    v_inner = compute_velocity(pos, dt)  # shape (N-1, 2)
    # Pad last row by repeating the final velocity
    return np.vstack([v_inner, v_inner[-1:]])


# =========================================================
# Convert per-droplet positions → flat row list for CSV
# =========================================================

def positions_to_rows(
    positions_by_id: Dict[int, List[List[float]]],
    dt: float = 1.0,
) -> List[List[float]]:
    """
    For each droplet, compute velocity and heading angle, then flatten into
    a list of rows ``[droplet_id, x, y, vx, vy, theta]``.

    Parameters
    ----------
    positions_by_id : dict
        Mapping droplet_id → list of [x, y] positions in frame order.
    dt : float
        Time step between consecutive frames.

    Returns
    -------
    rows : list of list
        Flat list of CSV rows, one per (droplet, frame) pair.
    """
    rows: List[List[float]] = []
    for droplet_id, positions in positions_by_id.items():
        pos = np.array(positions)
        if pos.size == 0:
            continue
        v = _velocity_for_csv(pos, dt)
        theta = np.arctan2(v[:, 1], v[:, 0])
        for idx in range(pos.shape[0]):
            rows.append([
                int(droplet_id),
                float(pos[idx, 0]),
                float(pos[idx, 1]),
                float(v[idx, 0]),
                float(v[idx, 1]),
                float(theta[idx]),
            ])
    return rows


# =========================================================
# Run YOLO tracking on videos and generate CSVs
# =========================================================

def generate_csvs(
    model: YOLO,
    videos_dir: Path,
    video_exts: Set[str],
    video_paths: List[Path],
    output_dir: Path = Path("output_droplet_tracking"),
    save_annotated: bool = True,
) -> None:
    """
    For every video in ``video_paths``:
      1. Run YOLO ``.track()`` to detect and track droplets.
      2. Collect bounding-box centres per tracked ID.
      3. Optionally save an annotated video.
      4. Write a CSV of positions/velocities.

    Parameters
    ----------
    model : YOLO
        Loaded Ultralytics YOLO model.
    videos_dir : Path
        Directory containing source videos.
    video_exts : set of str
        Accepted video file extensions (e.g. ``{".mp4", ".avi"}``).
    video_paths : list of Path
        Sorted list of video file paths to process.
    output_dir : Path
        Root output directory; sub-folders ``csv_files/`` and
        ``annotated_videos/`` are created automatically.
    save_annotated : bool
        If True, saves YOLO-annotated videos to ``annotated_videos/``.
    """
    if not video_paths:
        print(f"No videos found in {videos_dir.resolve()}")
        return

    annotated_dir = output_dir / "annotated_videos"
    csv_dir = output_dir / "csv_files"
    csv_dir.mkdir(parents=True, exist_ok=True)

    if save_annotated:
        annotated_dir.mkdir(parents=True, exist_ok=True)

    for video_path in video_paths:
        print(f"Tracking: {video_path.name}")

        track_results = model.track(
            source=str(video_path),
            stream=True,
            show=False,
            save=save_annotated,
            project=str(output_dir.resolve()),
            name="temp",
            exist_ok=True,
        )

        save_dir = None
        positions_by_id: Dict[int, List[List[float]]] = {}

        for r in track_results:
            if save_dir is None:
                save_dir = Path(r.save_dir)

            boxes = r.boxes
            if boxes.id is None:
                continue

            ids = boxes.id.cpu().numpy()
            xyxy = boxes.xyxy.cpu().numpy()

            for obj_id, box in zip(ids, xyxy):
                x1, y1, x2, y2 = box
                xc = (x1 + x2) / 2
                yc = (y1 + y2) / 2
                positions_by_id.setdefault(int(obj_id), []).append([xc, yc])

        # Move annotated video
        if save_annotated and save_dir and save_dir.exists():
            saved = list(save_dir.glob(f"*{video_path.suffix}"))
            if saved:
                target = annotated_dir / f"annotated_{video_path.name}"
                saved[0].replace(target)
                print(f"Annotated video saved to: {target}")

        # Write CSV
        rows = positions_to_rows(positions_by_id)
        csv_path = csv_dir / f"{video_path.stem}_droplets.csv"

        with csv_path.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["droplet id", "x", "y", "vx", "vy", "theta"])
            writer.writerows(rows)

        print(f"CSV saved to: {csv_path}")

    # Clean up temporary YOLO dir
    temp_dir = output_dir / "temp"
    if temp_dir.exists():
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)
