import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Ensure the project root is on sys.path
_PROJECT_ROOT = Path(__file__).resolve().parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from analysis.trajectories import load_positions_with_lengths
from analysis.plot_style import set_publication_style, StyleTokens, apply_grid, apply_legend

def plot_single_trajectory(positions: np.ndarray, droplet_id: int, output_dir: Path, experiment_name: str=""):
    """Plot and save a single droplet trajectory."""
    set_publication_style()
    
    fig, ax = plt.subplots(figsize=(6, 5))
    
    # Plot trajectory and mark start/end points
    ax.plot(positions[:, 0], positions[:, 1], lw=1.0, color=StyleTokens.TRAJECTORY, alpha=0.85)
    ax.plot(*positions[0], "o", ms=7, color=StyleTokens.START, label="Start", zorder=5)
    ax.plot(*positions[-1], "s", ms=7, color=StyleTokens.END, label="End", zorder=5)
    
    title = f"Trajectory - Droplet {droplet_id}"
    if experiment_name:
        title = f"{experiment_name} - " + title
        
    ax.set_title(title)
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    ax.set_aspect("equal")
    
    apply_legend(ax)
    apply_grid(ax, alpha=0.25)
    
    fig.tight_layout()
    
    # Save the plot
    prefix = f"{experiment_name}_" if experiment_name else ""
    output_path = output_dir / f"{prefix}droplet_{droplet_id}_trajectory.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Plot individual trajectories for each droplet.")
    parser.add_argument("csv_file", type=str, help="Path to the tracking CSV file")
    parser.add_argument("--output_dir", type=str, default="trajectory_plots", help="Directory to save the plots")
    
    args = parser.parse_args()
    csv_path = Path(args.csv_file)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Loading data from {csv_path}...")
    _, _, _, positions_dict, _ = load_positions_with_lengths(csv_path)
    
    experiment_name = csv_path.stem
    
    print(f"Found {len(positions_dict)} droplets. Generating plots...")
    for droplet_id, positions in positions_dict.items():
        if positions.shape[0] < 2:
            continue
        plot_single_trajectory(positions, droplet_id, output_dir, experiment_name=experiment_name)
        print(f"Saved plot for droplet {droplet_id}")
        
    print(f"All plots saved to {output_dir}")

if __name__ == "__main__":
    main()
