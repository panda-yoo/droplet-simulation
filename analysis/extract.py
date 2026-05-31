"""
analysis.extract — COMPATIBILITY SHIM
======================================
Functionality has moved to ``analysis.trajectories``.
This file re-exports both functions so that existing imports continue to work.

Do NOT add new implementations here.
"""

from analysis.trajectories import load_positions, load_positions_with_lengths

__all__ = ["load_positions", "load_positions_with_lengths"]
