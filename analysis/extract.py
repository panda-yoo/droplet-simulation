"""
analysis.extract
=================
Load droplet trajectory data from CSV files.

Contents
--------
- load_positions : read a tracking CSV → dict of per-droplet (N,2) arrays
"""

import numpy as np
from numpy.typing import NDArray
from typing import Dict, Tuple
import csv


# =========================================================
# CSV → position dict
# =========================================================

def load_positions(
    csv_path: str
) -> Tuple[Dict[int, NDArray[np.floating]], int, int]:
    """
    Read a CSV produced by either the simulation or YOLO tracking
    and return:
      - positions_consistent : dict mapping each consistent droplet
        ID to an (N, 2) array of (x, y) positions
      - total_paths          : total number of unique droplet IDs found
      - consistent_paths     : number of IDs kept (full-length tracks)

    Only droplets whose trajectory length equals the maximum
    observed length are kept (filters out partial tracks).
    """

    rows = []

    with open(csv_path, 'r') as f:

        reader = csv.reader(f)

        next(reader)  # skip header

        for row in reader:

            rows.append(row)

    # ── group positions by droplet ID ────────────────────
    positions_by_id = {}

    for row in rows:

        droplet_id = int(row[0])

        x = float(row[1])

        y = float(row[2])

        positions_by_id.setdefault(
            droplet_id,
            []
        ).append([x, y])

    # ── convert lists → numpy arrays ─────────────────────
    positions_by_id = {
        k: np.array(v)
        for k, v in positions_by_id.items()
    }

    # ── total number of paths found ──────────────────────
    total_paths = len(positions_by_id)

    # ── keep only droplets with full-length tracks ───────
    max_len = max(
        (
            v.shape[0]
            for v in positions_by_id.values()
        ),
        default=0
    )

    consistent_ids = [
        k for k, v in positions_by_id.items()
        if v.shape[0] == max_len
    ]

    positions_consistent = {
        k: positions_by_id[k]
        for k in consistent_ids
    }

    consistent_paths = len(positions_consistent)

    return positions_consistent, total_paths, consistent_paths
