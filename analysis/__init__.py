# =========================================================
# analysis package
# ---------------------------------------------------------
# Canonical analysis pipeline for droplet trajectory data.
#
#   analysis.types        — shared type aliases
#   analysis.trajectories — CSV loading / position extraction
#   analysis.msd          — MSD, fit, local exponent, pair-weighted
#   analysis.vacf         — VACF, velocity, persistence fit
#   analysis.orientation  — orientation correlation, persistence fit
#   analysis.psd          — PSD, power-law fit
#   analysis.diagnostics  — ensemble bias diagnostics
#   analysis.plotting     — all plot generation
#   analysis.main         — orchestration (post_processing, run_post_processing)
# =========================================================

from analysis.trajectories import load_positions, load_positions_with_lengths
