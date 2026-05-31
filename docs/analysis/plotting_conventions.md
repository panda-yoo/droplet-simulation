# Plotting Conventions

This document outlines the standard plotting philosophy and conventions for the entire active droplet project, ensuring thesis-quality visualizations across experimental analysis and simulations.

## Plotting Philosophy

All visualizations must follow a single, unified style to ensure that experimental summaries, ensemble averages, and simulation results can be presented side-by-side without visual inconsistency.

The project uses a centralized style module: `analysis/plot_style.py`. 

### Key principles:
1. **No ad-hoc styling**: Individual scripts should not set `matplotlib.rcParams` or use hardcoded string colors (like `"red"` or `"blue"`) unless they import from `StyleTokens`.
2. **Consistent Typography**: Sans-serif fonts, 11pt default, 12pt titles, and 9pt legends/ticks.
3. **Deterministic Coloring**: A specific droplet (e.g., Droplet 2) must always be plotted with the exact same color across *all* plots (Trajectory, MSD, VACF, etc.), whether in experimental data or simulations.

## Color Conventions

### Ensemble Curves & Highlights
`analysis.plot_style.StyleTokens` defines standard hex colors for ensemble-level or single-curve metrics to match the polished single-droplet summaries:

* **Trajectory Line**: `#2563eb` (Blue)
* **Start Marker**: `#16a34a` (Green)
* **End Marker**: `#dc2626` (Red)
* **MSD / PSD X**: `#2563eb` (Blue)
* **PSD Y**: `#dc2626` (Red)
* **Local Alpha ($\alpha$)**: `#7c3aed` (Purple)
* **VACF**: `#0891b2` (Cyan)
* **Orientation Correlation**: `#d97706` (Orange)
* **Fits**: `k:` (Black dotted line)

### Droplet Identification
When plotting multiple droplets on the same axes, use `get_droplet_color(droplet_id)` from `plot_style.py`. This uses a deterministic hash mapping to the `tab10` categorical palette, ensuring visual consistency across all comparisons.

## Figure Conventions
* **Grids**: All plots should use a faint background grid: `ax.grid(True, alpha=0.25)`. A helper function `apply_grid(ax)` is provided.
* **Line Widths**: Default curve line width is `1.5`. Auxiliary/fit lines should match or be slightly thicker (`2.0`) with dashed/dotted styles.
* **Markers**: Trajectory plots should explicitly mark the start with a circle (`"o"`, ms=6 or 7) and the end with a square (`"s"`, ms=6 or 7).

## Simulation Conventions
Simulation plots must mimic experimental analysis figures exactly:
* Trajectory plots use the same grid settings and `get_droplet_color` mapping.
* Start and End markers are explicitly added.
* Axis bounds should be set to `(0, 1)` with `aspect="equal"` to represent the normalized simulation domain.
