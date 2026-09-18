"""Reproduce Figure 6: a lower-boundary-to-interior constrained Brownian bridge.

Run: python3 Figure_006.py
Dependencies: numpy, scipy, matplotlib, and the adjacent ICB_gh_samplers.py.
The same 10,000 joint draws supply both marginal and sideways histograms.
The five additional paths use sequential conditional bridge draws.
"""

import csv
import math
from pathlib import Path
import random
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from ICB_gh_samplers import ICB_g_series, ICB_h_series, ICB_interpolate_step


SEED = 2026091606


def main():
    T = 3.0
    a = 2.0
    x = 0.0
    z = 1.8
    observation_times = [0.2, 2.8]
    t1, t2 = observation_times
    B = 10000
    seed = SEED
    path_seed = SEED + 1000
    print(f"Figure 6: T={T}, a={a}, x={x}, z={z}, B={B}, seed={seed}")

    # Each pair is sampled jointly: its second draw conditions on the first.
    random.seed(seed)
    started = time.perf_counter()
    first_values = []
    second_values = []
    for draw in range(B):
        first = ICB_interpolate_step(t1, T, x, z, a)
        second = ICB_interpolate_step(t2 - t1, T - t1, first, z, a)
        first_values.append(first)
        second_values.append(second)
    print(f"Generated {B} joint pairs in {time.perf_counter() - started:.3f} s")

    # Five separate paths on the 301-point grid, including exact endpoints.
    path_times = []
    for i in range(301):
        path_times.append(T * i / 300.0)
    for t in observation_times:
        if t not in path_times:
            path_times.append(t)
    path_times.sort()
    random.seed(path_seed)
    paths = []
    for path_number in range(5):
        values = [x]
        current_time = 0.0
        current_value = x
        for next_time in path_times[1:-1]:
            next_value = ICB_interpolate_step(
                next_time - current_time, T - current_time, current_value, z, a
            )
            values.append(next_value)
            current_time = next_time
            current_value = next_value
        values.append(z)
        paths.append(values)
    print(f"Generated 5 paths at {len(path_times)} times; path seed={path_seed}")

    # The marginal target is h(y;t,a) * g(y;T-t,z,a).
    # Scaled series cancel time-only exponential factors in normalization.
    y_grid = []
    for i in range(401):
        y_grid.append(a * i / 400.0)
    theoretical = []
    for t in observation_times:
        shape = []
        for y in y_grid:
            if y == 0.0 or y == a:
                shape.append(0.0)
                continue
            left = ICB_h_series(y, t, a, scaled=True)
            right = ICB_g_series(y, T - t, z, a, scaled=True)
            value = left * right
            if not math.isfinite(value) or value < -1e-12:
                raise FloatingPointError("Invalid theoretical density value")
            shape.append(max(0.0, value))
        normalizer = 0.0
        for i in range(1, len(y_grid)):
            normalizer += (shape[i - 1] + shape[i]) * (y_grid[i] - y_grid[i - 1]) / 2.0
        if not math.isfinite(normalizer) or normalizer <= 0.0:
            raise FloatingPointError("Invalid theoretical density normalizer")
        density = []
        for value in shape:
            density.append(value / normalizer)
        theoretical.append(density)

    samples = [first_values, second_values]
    kde_values = []
    for values in samples:
        kde = gaussian_kde(values)
        kde_values.append(kde(y_grid))

    fig = plt.figure(figsize=(8, 6))
    fig.subplots_adjust(left=0.085, right=0.975, bottom=0.09,
                        top=0.85, hspace=0.53, wspace=0.28)
    fig.suptitle(rf"Boundary-to-interior bridge: $x={x:g}$, $z={z:g}$, $a={a:g}$, $T={T:g}$",
                 fontsize=12, y=0.975)
    for panel in range(2):
        ax = fig.add_subplot(2, 2, panel + 1)
        ax.hist(samples[panel], bins=30, range=(0.0, a), density=True,
                color="0.65", alpha=0.6, edgecolor="black", linewidth=0.3,
                label="Histogram")
        ax.plot(y_grid, kde_values[panel], color="blue", linewidth=1.3, label="KDE")
        ax.plot(y_grid, theoretical[panel], color="red", linewidth=1.3, label="Theory")
        ax.set(xlim=(0.0, a), ylim=(0.0, None), xlabel=r"$y$", ylabel="Density",
               title=rf"$t={observation_times[panel]:.1f}$")
        ax.grid(alpha=0.15)
        ax.set_axisbelow(True)

    handles, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.945),
               frameon=False, ncol=3, fontsize=9)

    ax = fig.add_subplot(2, 1, 2)
    colors = ["#0072B2", "#009E73", "#E69F00", "#CC79A7", "#6A3D9A"]
    for path_number in range(5):
        ax.plot(path_times, paths[path_number], color=colors[path_number],
                linewidth=0.8, alpha=0.9)
        ax.plot(observation_times, [paths[path_number][path_times.index(t1)],
                                    paths[path_number][path_times.index(t2)]],
                "o", color=colors[path_number], markersize=3, zorder=5)
    ax.axhline(0.0, color="0.4", linestyle="--", linewidth=0.7)
    ax.axhline(a, color="0.4", linestyle="--", linewidth=0.7)
    right_limit = T + 0.85
    for panel in range(2):
        t = observation_times[panel]
        heights, edges = np.histogram(samples[panel], bins=30, range=(0.0, a), density=True)
        centers = (edges[:-1] + edges[1:]) / 2.0
        ax.barh(centers, 0.5 * heights, height=np.diff(edges), left=t,
                color="0.75", alpha=0.45, edgecolor="0.4", linewidth=0.3,
                label="Histogram" if panel == 0 else None)
        ax.plot(t + 0.5 * kde_values[panel], y_grid, color="blue", linewidth=1.25,
                label="KDE" if panel == 0 else None)
        ax.plot(t + 0.5 * np.array(theoretical[panel]), y_grid, color="red", linewidth=1.25,
                label="Theory" if panel == 0 else None)
        ax.axvline(t, color="0.55", linestyle=":", linewidth=0.8)
        right_limit = max(right_limit, t + 0.5 * max(max(heights),
                          max(kde_values[panel]), max(theoretical[panel])) + 0.1)
    ax.plot([0.0, T], [x, z], "o", color="black", markersize=4.5, zorder=6)
    ax.set(xlim=(-0.15, right_limit), ylim=(-0.1, a + 0.1),
           xlabel="Time", ylabel="Bridge value",
           title="Five bridge paths and marginal densities")
    ax.set_xticks([0.0, t1, T / 3.0, 2.0 * T / 3.0, t2, T])
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(alpha=0.15)
    ax.set_axisbelow(True)

    output = Path(__file__).resolve().parent / "figures"
    data_output = output / "data"
    data_output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        path = output / f"Figure_006.{extension}"
        fig.savefig(path, dpi=300)
        print(f"Saved {path}")
    plt.close(fig)
    path = data_output / "Figure_006_samples.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["draw", "y_t1", "y_t2"])
        for draw in range(B):
            writer.writerow([draw + 1, first_values[draw], second_values[draw]])
    print(f"Saved {path}")
    path = data_output / "Figure_006_paths.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path", "time", "value"])
        for path_number in range(5):
            for i in range(len(path_times)):
                writer.writerow([path_number + 1, path_times[i], paths[path_number][i]])
    print(f"Saved {path}")

    path = data_output / "Figure_006_density.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["t", "y", "theory_density", "kde_density"])
        for panel in range(2):
            for i in range(len(y_grid)):
                writer.writerow([observation_times[panel], y_grid[i],
                                 theoretical[panel][i], kde_values[panel][i]])
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
