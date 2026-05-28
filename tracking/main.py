"""
tracking.main
=============
Entry point for YOLO-based droplet tracking from video files.

Call ``run_tracking()`` to:
  1. Load the YOLO model (best.pt)
  2. Find all videos in the videos directory
  3. Generate per-video tracking CSVs
  4. Run post-processing analysis and plotting
"""

from ultralytics import YOLO
from pathlib import Path
from typing import Optional

import os

from tracking.utils import generate_csvs
from analysis.main import run_post_processing


# =========================================================
# Main tracking entry point
# =========================================================

def run_tracking(
    output_dir: str = "output_droplet_tracking",
    save_annotated: bool = True
):
    """
    Full tracking pipeline:
      1. Load YOLO model
      2. Discover video files
      3. Run tracking → generate CSVs
      4. Post-process and plot results
    """

    # =====================================================
    # PATHS
    # =====================================================

    BASE_DIR = Path(__file__).resolve().parents[1]

    # Resolve output_dir relative to BASE_DIR if it's relative
    output_path = Path(output_dir)
    if not output_path.is_absolute():
        output_path = (BASE_DIR / output_path).resolve()

    videos_dir = BASE_DIR / "videos"
    weight = BASE_DIR / "weight" / "best.pt"
    # =====================================================
    # VIDEO EXTENSIONS
    # =====================================================

    video_exts: set[str] = {
        ".mp4",
        ".mov",
        ".avi",
        ".mkv"
    }

    # =====================================================
    # FPS
    # =====================================================

    fps: float = 30.0

    dt: float = 1.0 / fps

    # =====================================================
    # LOAD MODEL
    # =====================================================

    model: Optional[YOLO] = None
    print(weight)

    try:
        model = YOLO(weight)
    except FileNotFoundError:

        print("model best.pt not found")

        return

    # =====================================================
    # VIDEO PATHS
    # =====================================================

    video_paths: list[Path] = sorted(
        [
            p for p in videos_dir.iterdir()
            if p.suffix.lower() in video_exts
        ]
    )

    # =====================================================
    # GENERATE CSVs
    # =====================================================

    generate_csvs(
        model=model,
        videos_dir=videos_dir,
        video_exts=video_exts,
        video_paths=video_paths,
        output_dir=output_path,
        save_annotated=save_annotated
    )

    # =====================================================
    # POST PROCESSING
    # =====================================================

    csv_dir = output_path / "csv_files"
    plots_dir = output_path / "plots"

    run_post_processing(
        base_dir=str(csv_dir),
        output_dir=str(plots_dir),
        dt=dt
    )
