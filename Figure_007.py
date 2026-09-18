"""Generate Figure 7 using simple loops and the current bridge samplers.

Run: python3 Figure_007.py
Dependencies: numpy, scipy, matplotlib, and the adjacent ICB_gh_samplers.py.
Outputs: PNG/SVG plots and CSV files containing joint draws and five paths.
The sampler automatically chooses its proposal orientation by time reversal.
"""

import csv
import math
from pathlib import Path
import random

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from ICB_gh_samplers import ICB_g_series, ICB_h_series, ICB_interpolate_step


# Edit the experiment settings here. Endpoints must stay in this figure's case.
A = 2.0
T = 3.0
X = 0.5
Z = A
OBSERVATION_TIMES = [0.2, 2.8]
N_DRAWS = 10000
N_PATHS = 5
PATH_STEPS = 300
GRID_POINTS = 401
SERIES_TERMS = 200
BINS = 30
SEED = 7001
SIDEWAYS_SCALE = 0.5


def theoretical_density(t, grid):
    """Normalize g(y;t,x,a)*h(a-y;T-t,a) by the trapezoidal rule."""
    values = []
    for y in grid:
        if y == 0.0 or y == A:
            value = 0.0
        else:
            # Scaling each series by exp(alpha) cancels on normalization.
            left = ICB_g_series(y, t, X, A, Nmax=SERIES_TERMS, scaled=True)
            right = ICB_h_series(A - y, T - t, A,
                                 Nmax=SERIES_TERMS, scaled=True)
            value = left * right
        if not math.isfinite(value) or value < -1e-12:
            raise ArithmeticError("Invalid theoretical density on the plot grid")
        values.append(max(0.0, value))

    normalizer = 0.0
    for j in range(1, len(grid)):
        normalizer += 0.5 * (values[j - 1] + values[j]) * (grid[j] - grid[j - 1])
    if normalizer <= 0.0:
        raise ArithmeticError("The density normalizer must be positive")
    density = []
    for value in values:
        density.append(value / normalizer)
    return density


def main():
    random.seed(SEED)
    t1, t2 = OBSERVATION_TIMES
    first_draws = []
    second_draws = []
    print(f"Figure 7: a={A:g}, T={T:g}, x={X:g}, z={Z:g}; seed={SEED}", flush=True)

    # Sequential conditioning produces joint pairs, not independent marginals.
    for draw in range(N_DRAWS):
        y1 = ICB_interpolate_step(t1, T, X, Z, A)
        y2 = ICB_interpolate_step(t2 - t1, T - t1, y1, Z, A)
        if not (0.0 < y1 < A and 0.0 < y2 < A):
            raise ArithmeticError("A sampled bridge value is outside the interval")
        first_draws.append(y1)
        second_draws.append(y2)
    samples = [first_draws, second_draws]
    print(f"Generated {N_DRAWS} joint pairs at times {t1:g} and {t2:g}.", flush=True)

    grid = np.linspace(0.0, A, GRID_POINTS)
    densities = []
    estimates = []
    for panel in range(2):
        densities.append(theoretical_density(OBSERVATION_TIMES[panel], grid))
        estimates.append(gaussian_kde(samples[panel])(grid))

    # Every path value is conditional on the previous value and final endpoint.
    # Form times by division, not repeated addition; include the marked times.
    path_times = []
    for j in range(PATH_STEPS + 1):
        path_times.append(T * j / PATH_STEPS)
    for t in OBSERVATION_TIMES:
        if t not in path_times:
            path_times.append(t)
    path_times.sort()
    random.seed(SEED + 1000)
    paths = []
    for path_number in range(N_PATHS):
        values = [X]
        previous_time = 0.0
        previous_value = X
        for j in range(1, len(path_times) - 1):
            t = path_times[j]
            value = ICB_interpolate_step(t - previous_time, T - previous_time,
                                         previous_value, Z, A)
            if not 0.0 < value < A:
                raise ArithmeticError("A path value is outside the interval")
            values.append(value)
            previous_time = t
            previous_value = value
        values.append(Z)
        paths.append(values)

    fig = plt.figure(figsize=(8, 6))
    axes = [fig.add_subplot(2, 2, 1), fig.add_subplot(2, 2, 2)]
    ax_path = fig.add_subplot(2, 1, 2)
    for panel in range(2):
        ax = axes[panel]
        ax.hist(samples[panel], bins=BINS, range=(0.0, A), density=True,
                color="0.65", alpha=0.6, edgecolor="black", linewidth=0.3)
        ax.plot(grid, estimates[panel], color="blue", linewidth=1.2)
        ax.plot(grid, densities[panel], color="red", linewidth=1.2)
        ax.set(xlim=(0.0, A), xlabel="y", ylabel=r"$k_{(3)}(y)$",
               title=f"t = {OBSERVATION_TIMES[panel]:g}")
        ax.set_axisbelow(True)
        ax.grid(color="0.92", linewidth=0.6)

    # Sideways marginal histograms, KDEs and densities share the joint samples.
    right_limit = T + 0.85
    for panel in range(2):
        t = OBSERVATION_TIMES[panel]
        heights, edges = np.histogram(samples[panel], bins=BINS,
                                      range=(0.0, A), density=True)
        for j in range(BINS):
            label = "Histogram" if panel == 0 and j == 0 else None
            ax_path.barh(edges[j], SIDEWAYS_SCALE * heights[j],
                         height=edges[j + 1] - edges[j], left=t, align="edge",
                         color="0.75", edgecolor="0.4", linewidth=0.3,
                         alpha=0.45, label=label)
        kde_x = []
        theory_x = []
        for j in range(len(grid)):
            kde_x.append(t + SIDEWAYS_SCALE * estimates[panel][j])
            theory_x.append(t + SIDEWAYS_SCALE * densities[panel][j])
        ax_path.plot(kde_x, grid, color="blue", linewidth=1.2,
                     label="KDE" if panel == 0 else None)
        ax_path.plot(theory_x, grid, color="red", linewidth=1.2,
                     label="Theory" if panel == 0 else None)
        right_limit = max(right_limit, max(kde_x) + 0.1, max(theory_x) + 0.1,
                          t + SIDEWAYS_SCALE * max(heights) + 0.1)
        ax_path.axvline(t, color="0.5", linestyle="--", linewidth=0.7)

    colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]
    for path_number in range(N_PATHS):
        color = colors[path_number % len(colors)]
        ax_path.plot(path_times, paths[path_number], color=color,
                     linewidth=0.7, alpha=0.7)
        for t in OBSERVATION_TIMES:
            j = path_times.index(t)
            ax_path.plot(t, paths[path_number][j], "o", color=color, markersize=4)
    ax_path.plot([0.0, 0.0, T, T, 0.0], [0.0, A, A, 0.0, 0.0],
                 color="darkgreen", linewidth=1.0)
    ax_path.plot([0.0, T], [X, Z], "o", color="black", markersize=4)
    ax_path.set(xlim=(-0.15, right_limit), ylim=(-0.1, A + 0.1),
                xlabel="t", ylabel="W(t)",
                title="Five bridge paths and marginal densities")
    fig.suptitle(f"Interior-to-boundary bridge: x={X:g}, z={Z:g}, a={A:g}, T={T:g}",
                 fontsize=11)
    handles, labels = ax_path.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.95),
               ncol=3, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))

    output = Path(__file__).resolve().parent / "figures"
    data_dir = output / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        path = output / f"Figure_007.{extension}"
        fig.savefig(path, dpi=300)
        print(f"Saved {path}")
    plt.close(fig)
    with (data_dir / "Figure_007_density.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["t", "y", "theory_density", "kde_density"])
        for panel in range(2):
            for j in range(len(grid)):
                writer.writerow([OBSERVATION_TIMES[panel], grid[j],
                                 densities[panel][j], estimates[panel][j]])
    with (data_dir / "Figure_007_samples.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["draw", "y_t1", "y_t2"])
        for draw in range(N_DRAWS):
            writer.writerow([draw + 1, first_draws[draw], second_draws[draw]])
    with (data_dir / "Figure_007_paths.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["path", "time", "value"])
        for path_number in range(N_PATHS):
            for j in range(len(path_times)):
                writer.writerow([path_number + 1, path_times[j], paths[path_number][j]])
    print(f"Saved joint draws and {N_PATHS} paths in {data_dir}")


if __name__ == "__main__":
    main()
