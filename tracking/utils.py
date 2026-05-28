"""
tracking.utils
==============
Helper functions for YOLO-based droplet tracking.

Contents
--------
- cal_velocity      : central-difference velocity from position arrays
- positions_to_rows : convert per-droplet positions → flat CSV rows
- generate_csvs     : run YOLO tracking on videos and write CSVs
"""

import numpy as np
from numpy.typing import NDArray
from typing import Dict, List, Set

from ultralytics import YOLO
from pathlib import Path
import csv


# =========================================================
# Velocity computation (central difference)
# =========================================================

def cal_velocity(
    arr: NDArray[np.floating],
    dt: float = 1.0
) -> NDArray[np.floating]:
    """
    Compute velocity from a (N, 2) position array using
    central differences (interior), forward (first), and
    backward (last) finite differences.
    """
    x = arr[:, 0]
    y = arr[:, 1]
    v = np.zeros_like(arr)

    if arr.shape[0] < 2:
        return v

    # central difference (interior points)
    v[1:-1, 0] = (x[2:] - x[:-2]) / (2 * dt)
    v[1:-1, 1] = (y[2:] - y[:-2]) / (2 * dt)

    # forward difference (first point)
    v[0, 0] = (x[1] - x[0]) / dt
    v[0, 1] = (y[1] - y[0]) / dt

    # backward difference (last point)
    v[-1, 0] = (x[-1] - x[-2]) / dt
    v[-1, 1] = (y[-1] - y[-2]) / dt

    return v


# =========================================================
# Convert per-droplet positions → flat row list for CSV
# =========================================================

def positions_to_rows(
    positions_by_id: Dict[int, List[List[float]]]
) -> List[List[float]]:
    """
    For each droplet, compute velocity and heading angle,
    then flatten into a list of rows:
      [droplet_id, x, y, vx, vy, theta]
    """
    rows = []
    for droplet_id, positions in positions_by_id.items():
        pos = np.array(positions)
        if pos.size == 0:
            continue
        v = cal_velocity(pos)
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
    output_dir: Path = Path('output_droplet_tracking'),
    save_annotated: bool = True
) -> None:
    """
    For every video in *video_paths*:
      1. Run YOLO .track() to detect and track droplets
      2. Collect bounding-box centres per tracked ID
      3. Optionally save an annotated video, and always write a CSV of positions/velocities
    """

    if not video_paths:
        print(f'No videos found in {videos_dir.resolve()}')

    else:

        # Create common output folders
        annotated_dir = output_dir / 'annotated_videos'
        csv_dir = output_dir / 'csv_files'

        csv_dir.mkdir(parents=True, exist_ok=True)
        if save_annotated:
            annotated_dir.mkdir(parents=True, exist_ok=True)

        for video_path in video_paths:

            print(f'Tracking: {video_path.name}')

            track_results = model.track(
                source=str(video_path),
                stream=True,
                show=False,
                save=save_annotated,
                project=str(output_dir.resolve()),
                name='temp',
                exist_ok=True,
            )

            save_dir = None
            positions_by_id = {}

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

                    positions_by_id.setdefault(
                        int(obj_id), []
                    ).append([xc, yc])

            # Move annotated video
            if save_annotated and save_dir and save_dir.exists():

                saved = list(
                    save_dir.glob(f'*{video_path.suffix}')
                )

                if saved:

                    target = (
                        annotated_dir /
                        f'annotated_{video_path.name}'
                    )

                    saved[0].replace(target)

                    print(
                        f'Annotated video saved to: {target}'
                    )

            # Save CSV
            rows = positions_to_rows(positions_by_id)

            csv_path = (
                csv_dir /
                f'{video_path.stem}_droplets.csv'
            )

            with csv_path.open('w', newline='') as f:

                writer = csv.writer(f)

                writer.writerow(
                    ['droplet id', 'x', 'y', 'vx', 'vy', 'theta']
                )

                writer.writerows(rows)

            print(f'CSV saved to: {csv_path}')

        # Clean up temporary YOLO dir if it was created
        temp_dir = output_dir / 'temp'
        if temp_dir.exists():
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)
