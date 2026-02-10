#!/usr/bin/env python3
"""Generate detailed seed-drop simulation visualizations."""
import csv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
import numpy as np

def load_trajectory_csv(path):
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

def load_seeds_csv(path):
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

def draw_tube_outline(ax, tube_type, projection='xz'):
    """Draw tube wall outline on the plot."""
    if tube_type == 'straight':
        r = 0.030
        ax.plot([-r, -r], [0, 0.30], 'k-', linewidth=2, alpha=0.7)
        ax.plot([r, r], [0, 0.30], 'k-', linewidth=2, alpha=0.7)
        ax.plot([-r, r], [0, 0], 'k-', linewidth=2, alpha=0.7)
        ax.plot([-r, r], [0.30, 0.30], 'k-', linewidth=1, alpha=0.3, linestyle='--')
    elif tube_type == 'curved':
        # Bent cylinder: bend_radius=0.573, bend_angle=30deg, center=(0,0,0.30)
        br = 0.573
        r = 0.030
        angle = 0.5236  # 30 degrees
        arc_center_x = br
        arc_center_z = 0.30
        thetas = np.linspace(0, angle, 100)
        # Inner and outer wall of bend
        for sign in [-1, 1]:
            wall_x = []
            wall_z = []
            for th in thetas:
                cx = arc_center_x + -(br + sign * r) * np.cos(th)
                cz = arc_center_z + -(br + sign * r) * np.sin(th)
                wall_x.append(cx)
                wall_z.append(cz)
            ax.plot(wall_x, wall_z, 'k-', linewidth=2, alpha=0.7)
    elif tube_type == 'funnel':
        r_top = 0.035
        r_bot = 0.018
        ax.plot([-r_top, -r_bot], [0.30, 0], 'k-', linewidth=2, alpha=0.7)
        ax.plot([r_top, r_bot], [0.30, 0], 'k-', linewidth=2, alpha=0.7)
        ax.plot([-r_bot, r_bot], [0, 0], 'k-', linewidth=2, alpha=0.7)
        ax.plot([-r_top, r_top], [0.30, 0.30], 'k-', linewidth=1, alpha=0.3, linestyle='--')


def plot_seed_drop_animation(output_dir, tube_type, title, save_path):
    """Create a multi-frame visualization showing seeds dropping through the tube."""
    traj_path = os.path.join(output_dir, 'trajectories.csv')
    seeds_path = os.path.join(output_dir, 'seeds.csv')

    trajectories = load_trajectory_csv(traj_path)
    seeds_data = load_seeds_csv(seeds_path)

    if not trajectories:
        print(f"No trajectory data in {output_dir}")
        return

    # Find time range
    all_times = []
    for pts in trajectories.values():
        for p in pts:
            all_times.append(p[0])
    t_min, t_max = min(all_times), max(all_times)

    # Create figure with 6 time snapshots + full trajectory + stats
    fig = plt.figure(figsize=(20, 12))
    fig.suptitle(f'Seed Drop Simulation: {title}', fontsize=16, fontweight='bold', y=0.98)

    # Time snapshots
    snapshot_times = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25]

    for idx, snap_t in enumerate(snapshot_times):
        ax = fig.add_subplot(2, 4, idx + 1)
        draw_tube_outline(ax, tube_type)

        # Draw seed positions at this time
        for pid, pts in trajectories.items():
            # Find closest time point
            best = None
            for p in pts:
                if best is None or abs(p[0] - snap_t) < abs(best[0] - snap_t):
                    best = p
            if best is not None and abs(best[0] - snap_t) < 0.005:
                # Draw seed as circle
                circle = plt.Circle((best[1], best[3]), 0.005,
                                   color='orangered', alpha=0.7, zorder=5)
                ax.add_patch(circle)

            # Draw faded trail up to this time
            trail_x = [p[1] for p in pts if p[0] <= snap_t]
            trail_z = [p[3] for p in pts if p[0] <= snap_t]
            if len(trail_x) > 1:
                ax.plot(trail_x, trail_z, color='steelblue', alpha=0.15, linewidth=0.5)

        ax.set_xlim(-0.06, 0.12)
        ax.set_ylim(-0.02, 0.32)
        ax.set_aspect('equal')
        ax.set_title(f't = {snap_t:.2f} s', fontsize=10)
        ax.set_xlabel('X (m)', fontsize=8)
        if idx == 0:
            ax.set_ylabel('Z (m)', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.2)

    # Full trajectory panel (XZ)
    ax_full = fig.add_subplot(2, 4, 7)
    draw_tube_outline(ax_full, tube_type)

    cmap = plt.cm.plasma
    n_traj = len(trajectories)
    for i, (pid, pts) in enumerate(list(trajectories.items())[:50]):
        xs = [p[1] for p in pts]
        zs = [p[3] for p in pts]
        ts = [p[0] for p in pts]

        # Color by time
        if len(xs) > 1:
            points = np.array([xs, zs]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            t_norm = np.array(ts[:-1])
            if t_max > t_min:
                t_norm = (t_norm - t_min) / (t_max - t_min)
            lc = LineCollection(segments, cmap=cmap, alpha=0.6, linewidth=0.8)
            lc.set_array(t_norm)
            ax_full.add_collection(lc)

    ax_full.set_xlim(-0.06, 0.12)
    ax_full.set_ylim(-0.02, 0.32)
    ax_full.set_aspect('equal')
    ax_full.set_title('Full Trajectories (XZ, colored by time)', fontsize=10)
    ax_full.set_xlabel('X (m)', fontsize=8)
    ax_full.set_ylabel('Z (m)', fontsize=8)
    ax_full.tick_params(labelsize=7)
    ax_full.grid(True, alpha=0.2)

    # Stats panel
    ax_stats = fig.add_subplot(2, 4, 8)
    ax_stats.axis('off')

    n_seeds = len(seeds_data)
    exit_speeds = [d['exit_speed'] for d in seeds_data]
    lateral_vels = [d['exit_lateral_velocity'] for d in seeds_data]
    wall_hits = [d['wall_hits'] for d in seeds_data]
    times = [d['time_to_exit'] for d in seeds_data]

    stats_text = (
        f"Tube Design: {title}\n"
        f"{'─' * 32}\n"
        f"Seeds simulated:  {n_seeds}\n"
        f"\n"
        f"Exit speed:       {np.mean(exit_speeds):.2f} ± {np.std(exit_speeds):.2f} m/s\n"
        f"Lateral velocity: {np.mean(lateral_vels):.3f} ± {np.std(lateral_vels):.3f} m/s\n"
        f"Time to exit:     {np.mean(times):.3f} ± {np.std(times):.3f} s\n"
        f"Wall hits (mean): {np.mean(wall_hits):.1f}\n"
        f"\n"
        f"Expected free-fall speed:\n"
        f"  v = √(2gh) = √(2×9.81×0.3)\n"
        f"    = {np.sqrt(2*9.81*0.3):.2f} m/s\n"
        f"\n"
        f"Simulation vs theory:\n"
        f"  {np.mean(exit_speeds):.2f} / {np.sqrt(2*9.81*0.3):.2f}"
        f" = {np.mean(exit_speeds)/np.sqrt(2*9.81*0.3):.1%}"
    )

    ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes,
                 fontsize=10, verticalalignment='top', fontfamily='monospace',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow',
                          edgecolor='gray', alpha=0.9))
    ax_stats.set_title('Simulation Summary', fontsize=10)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


def plot_combined_comparison(save_path):
    """Create a side-by-side trajectory comparison of all three tubes."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 8))
    fig.suptitle('Seed Drop Simulation — Trajectory Comparison', fontsize=16, fontweight='bold')

    configs = [
        ('output/straight_tube', 'straight', 'Straight Tube (30mm)'),
        ('output/curved_tube', 'curved', 'Curved Tube (30° bend)'),
        ('output/funnel_tube', 'funnel', 'Funnel Tube (35→18mm)'),
    ]

    for ax, (out_dir, tube_type, title) in zip(axes, configs):
        traj_path = os.path.join(out_dir, 'trajectories.csv')
        if not os.path.exists(traj_path):
            continue
        trajectories = load_trajectory_csv(traj_path)

        draw_tube_outline(ax, tube_type)

        cmap = plt.cm.viridis
        n = len(trajectories)
        for i, (pid, pts) in enumerate(list(trajectories.items())[:50]):
            color = cmap(i / max(n, 1))
            xs = [p[1] for p in pts]
            zs = [p[3] for p in pts]
            ax.plot(xs, zs, color=color, alpha=0.4, linewidth=0.8)
            # Draw seed at final position
            if xs:
                ax.plot(xs[-1], zs[-1], 'o', color=color, markersize=3, alpha=0.6)

        ax.set_xlim(-0.06, 0.12)
        ax.set_ylim(-0.02, 0.32)
        ax.set_aspect('equal')
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Z (m)')
        ax.grid(True, alpha=0.2)

        # Add gravity arrow
        ax.annotate('', xy=(0.09, 0.05), xytext=(0.09, 0.12),
                    arrowprops=dict(arrowstyle='->', color='red', lw=2))
        ax.text(0.095, 0.085, 'g', color='red', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


if __name__ == '__main__':
    os.chdir('/home/user/FluidX3D-DEM/dem/build')

    # Individual detailed seed-drop visualizations
    plot_seed_drop_animation('output/straight_tube', 'straight',
                           'Straight Tube (30mm diameter)',
                           '/home/user/FluidX3D-DEM/dem/docs/images/seed_drop_straight.png')

    plot_seed_drop_animation('output/curved_tube', 'curved',
                           'Curved Tube (30° bend)',
                           '/home/user/FluidX3D-DEM/dem/docs/images/seed_drop_curved.png')

    plot_seed_drop_animation('output/funnel_tube', 'funnel',
                           'Funnel Tube (35mm → 18mm)',
                           '/home/user/FluidX3D-DEM/dem/docs/images/seed_drop_funnel.png')

    # Side-by-side comparison
    plot_combined_comparison('/home/user/FluidX3D-DEM/dem/docs/images/seed_drop_comparison.png')
