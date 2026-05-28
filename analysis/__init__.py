# =========================================================
# analysis package
# ---------------------------------------------------------
# Post-processing of droplet trajectory CSVs.
#   analysis.extract   — CSV loading / position extraction
#   analysis.compute   — MSD, VACF, PSD, speed, orientation
#   analysis.main      — plotting and orchestration
# =========================================================

from analysis.extract import load_positions
