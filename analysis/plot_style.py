"""
analysis.plot_style
===================
Centralized style configuration and color mapping for all project plotting.
Ensures consistency between single-droplet summaries, ensemble plotting, and simulation output.
"""

import matplotlib
import matplotlib.pyplot as plt

class StyleTokens:
    """Hex color codes matching the polished single-droplet analysis."""
    TRAJECTORY = "#2563eb"
    START = "#16a34a"
    END = "#dc2626"
    MSD = "#2563eb"
    FIT = "k:"
    FIT_COLOR = "k"
    LOCAL_ALPHA = "#7c3aed"
    VACF = "#0891b2"
    ORIENT_CORR = "#d97706"
    PSD_X = "#2563eb"
    PSD_Y = "#dc2626"

def set_publication_style():
    """Apply standard typography and formatting to matplotlib."""
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "lines.linewidth": 1.5,
        "figure.dpi": 150,
    })

def apply_grid(ax, alpha=0.25, which="both"):
    """Apply a standardized grid to a given matplotlib axis."""
    ax.grid(True, which=which, alpha=alpha)

def get_droplet_color(droplet_id) -> str:
    """
    Get a deterministic color for a given droplet ID.
    Uses modulo mapping to consistently pick from a color palette, 
    ensuring Droplet X is always the same color across all plots.
    """
    palette = plt.cm.tab10.colors
    if isinstance(droplet_id, int):
        idx = droplet_id % len(palette)
    else:
        idx = hash(str(droplet_id)) % len(palette)
    return palette[idx]

def apply_legend(ax, loc="best", bbox_to_anchor=None):
    """Apply standard legend styling."""
    if bbox_to_anchor:
        ax.legend(loc=loc, bbox_to_anchor=bbox_to_anchor, framealpha=0.9, edgecolor="0.8")
    else:
        ax.legend(loc=loc, framealpha=0.9, edgecolor="0.8")

def add_diffusive_line(ax, alpha=1.0, label="Diffusive (α=1)"):
    """Add a gray dashed line for diffusive scaling (α=1)."""
    ax.axhline(alpha, linestyle="--", color="gray", lw=1.5, label=label)

def add_ballistic_line(ax, alpha=2.0, label="Ballistic (α=2)"):
    """Add a gray dotted line for ballistic scaling (α=2)."""
    ax.axhline(alpha, linestyle=":", color="gray", lw=1.5, alpha=0.6, label=label)

def add_zero_line(ax):
    """Add a black dashed line at y=0."""
    ax.axhline(0, linestyle="--", color="k", lw=1.0)
