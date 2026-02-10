#!/usr/bin/env python3
"""
DEM Seed Drop Simulator - Postprocessing & Visualization
Reads CSV output from dem_sim and generates summary plots.

Usage:
    python dem_postprocess.py <output_dir> [--compare <dir2> [<dir3> ...]]
    python dem_postprocess.py output/straight_tube --compare output/curved_tube output/spiral_tube
"""

import argparse
import csv
import os
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not available. Install with: pip install matplotlib")

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def load_seeds_csv(path):
    """Load per-seed metrics CSV into a list of dicts."""
    data = []
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            entry = {}
            for k, v in row.items():
                try:
                    entry[k] = float(v)
                except ValueError:
                    entry[k] = v
            data.append(entry)
    return data


def load_summary_csv(path):
    """Load summary CSV into a dict."""
    summary = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                summary[row['metric']] = float(row['value'])
            except (ValueError, KeyError):
                pass
    return summary


def load_trajectory_csv(path):
    """Load trajectory CSV into dict of {particle_id: [(t, x, y, z), ...]}."""
    trajectories = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pid = int(float(row['particle_id']))
            t = float(row['time'])
            x = float(row['x'])
            y = float(row['y'])
            z = float(row['z'])
            if pid not in trajectories:
                trajectories[pid] = []
            trajectories[pid].append((t, x, y, z))
    return trajectories


def print_summary(name, summary):
    """Print summary to console."""
    print(f"\n{'='*55}")
    print(f"  {name}")
    print(f"{'='*55}")
    for k, v in summary.items():
        print(f"  {k:30s}: {v:.6f}")
    print(f"{'='*55}")


def plot_single_run(output_dir, save_path=None):
    """Generate plots for a single simulation run."""
    if not HAS_MATPLOTLIB:
        print("Cannot generate plots without matplotlib")
        return

    seeds_path = os.path.join(output_dir, 'seeds.csv')
    summary_path = os.path.join(output_dir, 'summary.csv')

    if not os.path.exists(seeds_path):
        print(f"Error: {seeds_path} not found")
        return

    data = load_seeds_csv(seeds_path)
    if not data:
        print("No seed data found")
        return

    summary = load_summary_csv(summary_path) if os.path.exists(summary_path) else {}

    # Extract arrays
    times = [d['time_to_exit'] for d in data]
    wall_hits = [d['wall_hits'] for d in data]
    exit_speeds = [d['exit_speed'] for d in data]
    lateral_vels = [d['exit_lateral_velocity'] for d in data]
    impulses = [d['cumulative_impulse'] for d in data]
    bounce_energy = [d['max_bounce_energy'] for d in data]

    run_name = os.path.basename(output_dir)

    fig = plt.figure(figsize=(16, 12))
    fig.suptitle(f'DEM Seed Drop Simulation: {run_name}', fontsize=14, fontweight='bold')
    gs = gridspec.GridSpec(3, 3, hspace=0.35, wspace=0.3)

    # 1. Time to exit histogram
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(times, bins=30, color='steelblue', edgecolor='black', alpha=0.8)
    ax1.set_xlabel('Time to exit (s)')
    ax1.set_ylabel('Count')
    ax1.set_title('Time to Exit Distribution')
    if 'mean_time_to_exit' in summary:
        ax1.axvline(summary['mean_time_to_exit'], color='red', linestyle='--',
                     label=f"mean={summary['mean_time_to_exit']:.4f}")
        ax1.legend(fontsize=8)

    # 2. Wall hits histogram
    ax2 = fig.add_subplot(gs[0, 1])
    max_hits = max(wall_hits) if wall_hits else 0
    bins = range(0, int(max_hits) + 2)
    ax2.hist(wall_hits, bins=bins, color='coral', edgecolor='black', alpha=0.8)
    ax2.set_xlabel('Wall hits')
    ax2.set_ylabel('Count')
    ax2.set_title('Wall Hits Distribution')

    # 3. Exit speed histogram
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.hist(exit_speeds, bins=30, color='seagreen', edgecolor='black', alpha=0.8)
    ax3.set_xlabel('Exit speed (m/s)')
    ax3.set_ylabel('Count')
    ax3.set_title('Exit Speed Distribution')

    # 4. Lateral velocity histogram
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.hist(lateral_vels, bins=30, color='mediumpurple', edgecolor='black', alpha=0.8)
    ax4.set_xlabel('Exit lateral velocity (m/s)')
    ax4.set_ylabel('Count')
    ax4.set_title('Exit Lateral Velocity Distribution')

    # 5. Cumulative impulse histogram
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.hist(impulses, bins=30, color='goldenrod', edgecolor='black', alpha=0.8)
    ax5.set_xlabel('Cumulative wall impulse (N·s)')
    ax5.set_ylabel('Count')
    ax5.set_title('Wall Impulse Distribution')

    # 6. Max bounce energy histogram
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.hist(bounce_energy, bins=30, color='indianred', edgecolor='black', alpha=0.8)
    ax6.set_xlabel('Max bounce energy (J)')
    ax6.set_ylabel('Count')
    ax6.set_title('Max Bounce Energy Distribution')

    # 7. Exit velocity scatter (vx vs vy)
    ax7 = fig.add_subplot(gs[2, 0])
    vx = [d['exit_vx'] for d in data]
    vy = [d['exit_vy'] for d in data]
    ax7.scatter(vx, vy, s=5, alpha=0.5, c='steelblue')
    ax7.set_xlabel('Exit Vx (m/s)')
    ax7.set_ylabel('Exit Vy (m/s)')
    ax7.set_title('Exit Velocity Scatter (XY)')
    ax7.set_aspect('equal')
    ax7.axhline(0, color='gray', linewidth=0.5)
    ax7.axvline(0, color='gray', linewidth=0.5)

    # 8. Time to exit vs wall hits
    ax8 = fig.add_subplot(gs[2, 1])
    ax8.scatter(times, wall_hits, s=5, alpha=0.5, c='coral')
    ax8.set_xlabel('Time to exit (s)')
    ax8.set_ylabel('Wall hits')
    ax8.set_title('Time to Exit vs Wall Hits')

    # 9. Summary text
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.axis('off')
    summary_text = f"Seeds: {len(data)}\n"
    if summary:
        summary_text += f"Time to exit: {summary.get('mean_time_to_exit', 0):.4f} ± {summary.get('sd_time_to_exit', 0):.4f} s\n"
        summary_text += f"Wall hits: {summary.get('mean_wall_hits', 0):.1f} ± {summary.get('sd_wall_hits', 0):.1f}\n"
        summary_text += f"Exit speed: {summary.get('mean_exit_speed', 0):.3f} ± {summary.get('sd_exit_speed', 0):.3f} m/s\n"
        summary_text += f"Lateral vel: {summary.get('mean_exit_lateral_vel', 0):.4f} ± {summary.get('sd_exit_lateral_vel', 0):.4f} m/s\n"
        summary_text += f"\nSpacing Risk Proxy: {summary.get('spacing_risk_proxy', 0):.4f}\n"
        summary_text += "(lower = better singulation)"
    ax9.text(0.1, 0.9, summary_text, transform=ax9.transAxes,
             fontsize=10, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    ax9.set_title('Summary')

    if save_path is None:
        save_path = os.path.join(output_dir, 'plots.png')
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {save_path}")
    plt.close()


def plot_comparison(dirs, labels=None, save_path=None):
    """Generate comparison plots for multiple simulation runs."""
    if not HAS_MATPLOTLIB:
        print("Cannot generate plots without matplotlib")
        return

    if labels is None:
        labels = [os.path.basename(d) for d in dirs]

    all_data = []
    all_summaries = []
    for d in dirs:
        seeds_path = os.path.join(d, 'seeds.csv')
        summary_path = os.path.join(d, 'summary.csv')
        if os.path.exists(seeds_path):
            all_data.append(load_seeds_csv(seeds_path))
        else:
            all_data.append([])
        if os.path.exists(summary_path):
            all_summaries.append(load_summary_csv(summary_path))
        else:
            all_summaries.append({})

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('DEM Seed Drop Simulation Comparison', fontsize=14, fontweight='bold')

    colors = ['steelblue', 'coral', 'seagreen', 'mediumpurple', 'goldenrod']

    metrics_to_plot = [
        ('time_to_exit', 'Time to Exit (s)', axes[0, 0]),
        ('wall_hits', 'Wall Hits', axes[0, 1]),
        ('exit_speed', 'Exit Speed (m/s)', axes[0, 2]),
        ('exit_lateral_velocity', 'Exit Lateral Velocity (m/s)', axes[1, 0]),
        ('cumulative_impulse', 'Wall Impulse (N·s)', axes[1, 1]),
    ]

    for metric, title, ax in metrics_to_plot:
        for i, (data, label) in enumerate(zip(all_data, labels)):
            if not data:
                continue
            values = [d[metric] for d in data]
            c = colors[i % len(colors)]
            ax.hist(values, bins=25, alpha=0.5, label=label, color=c, edgecolor='black', linewidth=0.5)
        ax.set_title(title)
        ax.set_ylabel('Count')
        ax.legend(fontsize=8)

    # Summary comparison bar chart
    ax_bar = axes[1, 2]
    bar_metrics = ['spacing_risk_proxy']
    x_pos = range(len(dirs))
    values = [s.get('spacing_risk_proxy', 0) for s in all_summaries]
    bar_colors = [colors[i % len(colors)] for i in range(len(dirs))]
    ax_bar.bar(x_pos, values, color=bar_colors, edgecolor='black')
    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    ax_bar.set_title('Spacing Risk Proxy (lower = better)')
    ax_bar.set_ylabel('Risk Proxy')

    plt.tight_layout()

    if save_path is None:
        save_path = 'comparison.png'
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved comparison plot: {save_path}")
    plt.close()


def plot_trajectories(output_dir, save_path=None):
    """Plot 3D trajectory visualization."""
    if not HAS_MATPLOTLIB:
        return

    traj_path = os.path.join(output_dir, 'trajectories.csv')
    if not os.path.exists(traj_path):
        print(f"No trajectory file found: {traj_path}")
        return

    trajectories = load_trajectory_csv(traj_path)
    if not trajectories:
        return

    fig = plt.figure(figsize=(12, 8))

    # 2D projections (XZ and YZ)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    cmap = plt.cm.viridis
    n_traj = len(trajectories)

    for i, (pid, points) in enumerate(list(trajectories.items())[:50]):
        color = cmap(i / max(n_traj, 1))
        xs = [p[1] for p in points]
        zs = [p[3] for p in points]
        ys = [p[2] for p in points]

        ax1.plot(xs, zs, color=color, alpha=0.5, linewidth=0.5)
        ax2.plot(ys, zs, color=color, alpha=0.5, linewidth=0.5)

    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Z (m)')
    ax1.set_title('Trajectories (XZ projection)')
    ax1.set_aspect('equal')

    ax2.set_xlabel('Y (m)')
    ax2.set_ylabel('Z (m)')
    ax2.set_title('Trajectories (YZ projection)')
    ax2.set_aspect('equal')

    run_name = os.path.basename(output_dir)
    fig.suptitle(f'Seed Trajectories: {run_name}', fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_path is None:
        save_path = os.path.join(output_dir, 'trajectories.png')
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved trajectory plot: {save_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='DEM Seed Drop Simulator - Postprocessing')
    parser.add_argument('output_dir', help='Path to simulation output directory')
    parser.add_argument('--compare', nargs='+', help='Additional output dirs for comparison')
    parser.add_argument('--no-plots', action='store_true', help='Skip plot generation')
    parser.add_argument('--save-path', help='Custom save path for plots')
    args = parser.parse_args()

    # Print summary
    summary_path = os.path.join(args.output_dir, 'summary.csv')
    if os.path.exists(summary_path):
        summary = load_summary_csv(summary_path)
        print_summary(os.path.basename(args.output_dir), summary)

    if args.compare:
        for d in args.compare:
            sp = os.path.join(d, 'summary.csv')
            if os.path.exists(sp):
                print_summary(os.path.basename(d), load_summary_csv(sp))

    if not args.no_plots and HAS_MATPLOTLIB:
        # Single run plots
        plot_single_run(args.output_dir, args.save_path)
        plot_trajectories(args.output_dir)

        # Comparison
        if args.compare:
            all_dirs = [args.output_dir] + args.compare
            plot_comparison(all_dirs)


if __name__ == '__main__':
    main()
