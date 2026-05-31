"""
analysis.trajectories
=====================
Load droplet trajectory data from CSV files.

Contents
--------
- load_positions              : read CSV → consistent (full-length) tracks only
- load_positions_with_lengths : read CSV → consistent + all tracks + lengths
"""

import csv
import numpy as np
from numpy.typing import NDArray
from pathlib import Path
from typing import Dict, Tuple

from analysis.types import PositionsDict


# =========================================================
# CSV → position dict
# =========================================================

def load_positions(
    csv_path: Path,
) -> Tuple[PositionsDict, int, int]:
    """
    Read a tracking CSV and return only the full-length (consistent) tracks.

    Parameters
    ----------
    csv_path : Path
        Path to a CSV with columns: droplet_id, x, y, [vx, vy, theta].

    Returns
    -------
    positions_consistent : PositionsDict
        Mapping droplet_id → (N, 2) array of (x, y) positions.
        Only droplets whose trajectory length equals the global maximum are kept.
    total_paths : int
        Total number of unique droplet IDs found in the CSV.
    consistent_paths : int
        Number of droplet IDs retained (full-length tracks only).

    Notes
    -----
    Partial tracks are discarded so that every trajectory in the returned dict
    has identical length, which is required for unbiased ensemble averaging.
    """
    rows: list[list[str]] = []
    with open(csv_path, "r") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            rows.append(row)

    positions_by_id: Dict[int, list[list[float]]] = {}
    for row in rows:
        droplet_id = int(row[0])
        x = float(row[1])
        y = float(row[2])
        positions_by_id.setdefault(droplet_id, []).append([x, y])

    positions_by_id_np: PositionsDict = {
        k: np.array(v) for k, v in positions_by_id.items()
    }

    total_paths = len(positions_by_id_np)

    max_len = max(
        (v.shape[0] for v in positions_by_id_np.values()),
        default=0,
    )

    consistent_ids = [
        k for k, v in positions_by_id_np.items() if v.shape[0] == max_len
    ]

    positions_consistent: PositionsDict = {
        k: positions_by_id_np[k] for k in consistent_ids
    }

    consistent_paths = len(positions_consistent)

    return positions_consistent, total_paths, consistent_paths


def load_positions_with_lengths(
    csv_path: Path,
) -> Tuple[PositionsDict, int, int, PositionsDict, Dict[int, int]]:
    """
    Read a CSV and return consistent tracks together with all tracks and lengths.

    Parameters
    ----------
    csv_path : Path
        Path to a CSV with columns: droplet_id, x, y, [vx, vy, theta].

    Returns
    -------
    positions_consistent : PositionsDict
        Full-length trajectories only.
    total_paths : int
        Total number of unique droplet IDs.
    consistent_paths : int
        Number of full-length trajectories.
    positions_all : PositionsDict
        All trajectories including partial tracks.
    lengths : dict
        Track length (number of frames) per droplet ID.

    Notes
    -----
    ``positions_consistent`` is a strict subset of ``positions_all``.
    ``lengths`` covers every droplet ID in ``positions_all``.
    """
    positions_consistent, total_paths, consistent_paths = load_positions(csv_path)

    rows: list[list[str]] = []
    with open(csv_path, "r") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            rows.append(row)

    positions_by_id: Dict[int, list[list[float]]] = {}
    for row in rows:
        droplet_id = int(row[0])
        x = float(row[1])
        y = float(row[2])
        positions_by_id.setdefault(droplet_id, []).append([x, y])

    positions_all: PositionsDict = {
        k: np.array(v) for k, v in positions_by_id.items()
    }

    lengths: Dict[int, int] = {k: v.shape[0] for k, v in positions_all.items()}

    return (
        positions_consistent,
        total_paths,
        consistent_paths,
        positions_all,
        lengths,
    )
