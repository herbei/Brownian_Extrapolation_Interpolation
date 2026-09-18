"""Figure 1: interior-start extrapolation in both Devroye sampling regimes.

Run this file from any working directory. Outputs are saved beside this
script in figures/, with the simulated values in figures/data/.
"""

import csv
import math
import random
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from ICB_gh_samplers import ICB_g_series, ICB_sample_g


# Edit these constants to change the experiment.
SEED = 1001
B = 10000
X = 0.5
A_VALUES = [1.8, 2.5]
T_VALUES = [1.0, 0.2]
GRID_N = 1001
N_TERMS = 500
HISTOGRAM_BINS = 30
DPI = 300


def main():
    random.seed(SEED)
    output_dir = Path(__file__).resolve().parent / "figures"
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(1, len(A_VALUES), figsize=(8, 3), squeeze=False)
    sample_rows = []
    density_rows = []

    for panel in range(len(A_VALUES)):
        a = A_VALUES[panel]
        t = T_VALUES[panel]
        samples = []
        for draw in range(B):
            value = ICB_sample_g(t, X, a)
            if not math.isfinite(value) or not 0.0 < value < a:
                raise ValueError("A simulated value is not finite and inside (0, a).")
            samples.append(value)
            sample_rows.append([panel + 1, draw + 1, t, X, a, value])

        grid = np.linspace(0.0, a, GRID_N)
        raw_density = []
        for y in grid:
            if y == 0.0 or y == a:
                value = 0.0
            else:
                # The common exp(-alpha) factor cancels during normalization.
                value = ICB_g_series(y, t, X, a, Nmax=N_TERMS, scaled=True)
            if not math.isfinite(value) or value < 0.0:
                raise ValueError("The theoretical density is not finite and nonnegative.")
            raw_density.append(value)

        normalizer = 0.0
        for j in range(1, GRID_N):
            width = grid[j] - grid[j - 1]
            normalizer += 0.5 * (raw_density[j - 1] + raw_density[j]) * width
        if not math.isfinite(normalizer) or normalizer <= 0.0:
            raise ValueError("The theoretical density has no positive finite integral.")

        density = []
        for value in raw_density:
            density.append(value / normalizer)

        integral = 0.0
        for j in range(1, GRID_N):
            width = grid[j] - grid[j - 1]
            integral += 0.5 * (density[j - 1] + density[j]) * width
        if abs(integral - 1.0) > 1e-10:
            raise ValueError("The normalized density does not integrate to one.")

        kde = gaussian_kde(samples)
        kde_density = kde(grid)
        for j in range(GRID_N):
            density_rows.append([panel + 1, t, X, a, grid[j], density[j], kde_density[j]])

        ax = axes[0, panel]
        ax.hist(samples, bins=HISTOGRAM_BINS, range=(0.0, a), density=True,
                color="0.65", edgecolor="black", linewidth=0.3, alpha=0.6)
        ax.plot(grid, kde_density, color="blue", linewidth=1.2, label="KDE")
        ax.plot(grid, density, color="red", linewidth=1.2, label="Theory")
        ratio = a / math.sqrt(t)
        if ratio < 2.0:
            regime = "< 2"
        elif ratio > 2.0:
            regime = "> 2"
        else:
            # Equality is already shown by the numerical ratio in the title.
            regime = ""
        ax.set_title(f"a = {a:g}, t = {t:g}, x = {X:g}\n"
                     + rf"$a/\sqrt{{t}} = {ratio:.3g} {regime}$", fontsize=10)
        ax.set_xlabel("y")
        ax.set_ylabel("Density")
        ax.set_xlim(0.0, a)
        ax.set_ylim(bottom=0.0)
        ax.grid(alpha=0.18)
        ax.set_axisbelow(True)
        ax.legend(fontsize=8, frameon=False)
        print(f"Figure 1: a={a:g}, t={t:g}, x={X:g}, a/sqrt(t)={ratio:g}, "
              f"samples={len(samples)}, density integral={integral:.12f}")

    with (data_dir / "Figure_001_samples.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["panel", "draw", "t", "x", "a", "y"])
        for row in sample_rows:
            writer.writerow(row)
    with (data_dir / "Figure_001_density.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["panel", "t", "x", "a", "y", "theory_density", "kde_density"])
        for row in density_rows:
            writer.writerow(row)

    figure.tight_layout()
    figure.savefig(output_dir / "Figure_001.png", dpi=DPI)
    figure.savefig(output_dir / "Figure_001.svg")
    plt.close(figure)
    print(f"Saved Figure_001.png, Figure_001.svg, and CSV data in {output_dir}")


if __name__ == "__main__":
    main()
