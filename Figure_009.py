"""Reproduce Figure 9: the two continuous truncation-bound curves.

Run: python3 Figure_009.py
Dependency: matplotlib (the calculations use only Python's math module).

This reproduces the continuous square-root expressions in manuscript Figure 9,
before the ceiling and minus-one operations in the draft's inequalities.
For integer N, the displayed condition N > ceil(q) - 1 means N >= ceil(q).
These expressions do not enforce the proposition's separate N0 condition.
The reversed beta axis is a display choice, not beta = pi**2/(4*alpha).
"""

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


TOL = 1e-16
PARAMETER_MIN = 0.001
PARAMETER_MAX = 10.0
GRID_POINTS = 500


def main():
    parameters = []
    alpha_thresholds = []
    beta_thresholds = []
    alpha_log10_n = []
    beta_log10_n = []
    beta_plot_y = []

    # Use the linear parameter grid and formulas in the manuscript figure.
    for i in range(GRID_POINTS):
        parameter = PARAMETER_MIN + (PARAMETER_MAX - PARAMETER_MIN) * i / (GRID_POINTS - 1)
        alpha_n = math.sqrt((2.0 / parameter) * math.log(4.0 / (parameter * TOL)))
        beta_n = math.sqrt((1.0 / parameter) * math.log(2.0 / TOL))
        parameters.append(parameter)
        alpha_thresholds.append(alpha_n)
        beta_thresholds.append(beta_n)
        alpha_log10_n.append(math.log10(alpha_n))
        beta_log10_n.append(math.log10(beta_n))
        beta_plot_y.append(PARAMETER_MIN + PARAMETER_MAX - parameter)

    fig, ax = plt.subplots(figsize=(6, 2.5))
    fig.subplots_adjust(left=0.10, right=0.62, bottom=0.22, top=0.96)
    ax.plot(alpha_log10_n, parameters, color="black", linewidth=1.3,
            label=r"$\log_{10}(N(\alpha))$")
    ax.plot(beta_log10_n, beta_plot_y, color="0.4", linestyle="--",
            linewidth=1.3, label=r"$\log_{10}(N(\beta))$")
    ax.set_xlabel(r"$\log_{10}(N)$")
    ax.set_ylabel(r"$\alpha$")
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax.yaxis.set_major_locator(MultipleLocator(2.5))
    ax.yaxis.set_minor_locator(MultipleLocator(1.25))
    ax.grid(which="major", color="0.90", linewidth=0.8)
    ax.grid(which="minor", color="0.94", linewidth=0.5)
    ax.tick_params(which="minor", length=0)

    # The two parameters range independently over [0.001, 10]. Reflect the
    # right axis so beta decreases upwards, exactly as in the original plot.
    ax_beta = ax.twinx()
    lower, upper = ax.get_ylim()
    axis_sum = PARAMETER_MIN + PARAMETER_MAX
    ax_beta.set_ylim(axis_sum - lower, axis_sum - upper)
    ax_beta.set_ylabel(r"$\beta$")
    ax_beta.yaxis.set_major_locator(MultipleLocator(2.5))
    ax.legend(loc="center left", bbox_to_anchor=(1.26, 0.5),
              frameon=False, fontsize=9)

    output = Path(__file__).resolve().parent / "figures"
    data_output = output / "data"
    data_output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        path = output / f"Figure_009.{extension}"
        fig.savefig(path, dpi=300)
        print(f"Saved {path}")
    plt.close(fig)

    path = data_output / "Figure_009_bounds.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["parameter", "alpha_continuous_N", "beta_continuous_N",
                         "alpha_log10_N", "beta_log10_N", "beta_plot_y",
                         "alpha_ceil_N", "beta_ceil_N"])
        for i in range(GRID_POINTS):
            writer.writerow([parameters[i], alpha_thresholds[i], beta_thresholds[i],
                             alpha_log10_n[i], beta_log10_n[i], beta_plot_y[i],
                             math.ceil(alpha_thresholds[i]), math.ceil(beta_thresholds[i])])
    print(f"Saved {path}")
    print(f"tol={TOL:g}; {GRID_POINTS} points from {PARAMETER_MIN:g} to {PARAMETER_MAX:g}")
    print(f"Continuous alpha thresholds: {min(alpha_thresholds):.6f} to {max(alpha_thresholds):.6f}")
    print(f"Continuous beta thresholds:  {min(beta_thresholds):.6f} to {max(beta_thresholds):.6f}")


if __name__ == "__main__":
    main()
