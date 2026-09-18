"""Reproduce Figure 3: the acceptance-rate curves A1(beta) and A2(beta).

Run: python3 Figure_003.py
Dependencies: numpy, matplotlib, and the adjacent ICB_gh_samplers.py.

Uses the verified cutoff constants shared with the samplers. The 200-term
sums are sufficient on the plotted range 0.1 <= beta <= 10. This plots the
paper's theoretical A1 and A2, not empirical sampler acceptance rates.
"""

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

from ICB_gh_samplers import ICB_ALPHA0, ICB_BETA0


def acceptance_a1(beta, terms=200):
    """Evaluate the Section 3.1 expression for A1 on beta >= beta0."""
    numerator = math.fsum([1.0] + [
        2.0 * (-1)**n * math.exp(-beta * n**2)
        for n in range(1, terms + 1)
    ])
    series = math.fsum(math.exp(-beta * n**2)
                       for n in range(1, terms + 1))
    denominator = -math.expm1(-beta) + beta * series
    return numerator / denominator


def acceptance_a2(beta, terms=200):
    """Evaluate A2 on 0 < beta <= beta0 using the positive-series ratio.

    The appendix gives A2 = (8/pi**2) * (sum over positive odd n of
    exp(-alpha*n**2)) / (sum over positive n of n**2*exp(-alpha*n**2)).
    Cancel exp(-alpha) to avoid tiny numerator and denominator values.
    """
    alpha = math.pi**2 / (4.0 * beta)
    numerator = math.fsum(math.exp(-alpha * (n**2 - 1))
                          for n in range(1, terms + 1, 2))
    denominator = math.fsum(n**2 * math.exp(-alpha * (n**2 - 1))
                            for n in range(1, terms + 1))
    return (8.0 / math.pi**2) * numerator / denominator


def main():
    beta0 = ICB_BETA0
    beta_left = np.linspace(0.1, beta0, 300)
    beta_right = np.linspace(beta0, 10.0, 300)
    a2 = [acceptance_a2(beta) for beta in beta_left]
    a1 = [acceptance_a1(beta) for beta in beta_right]

    print(f"alpha_0 = {ICB_ALPHA0:.12f}")
    print(f"beta_0  = {beta0:.12f}")
    print(f"A1(beta_0) = {acceptance_a1(beta0):.12f}")
    print(f"A2(beta_0) = {acceptance_a2(beta0):.12f}")

    fig, ax = plt.subplots(figsize=(6, 3))
    fig.subplots_adjust(left=0.12, right=0.71, bottom=0.20, top=0.96)
    ax.plot(beta_left, a2, color="red", linewidth=1.3,
            label=r"$A_2(\beta_t)$")
    ax.plot(beta_right, a1, color="blue", linewidth=1.3,
            label=r"$A_1(\beta_t)$")
    ax.axvline(beta0, color="0.7", linestyle="--", linewidth=1.0)
    ax.plot(beta0, acceptance_a2(beta0), "o", color="red", markersize=5)
    ax.text(beta0, 0.605, rf"$\beta_0={beta0:.6f}$", ha="center")
    ax.text(1.0, 0.69, r"$\beta_t<\beta_0$", color="red",
            ha="center", fontsize=12)
    ax.text(5.0, 0.69, r"$\beta_t>\beta_0$", color="blue",
            ha="center", fontsize=12)
    ax.set(xlim=(-0.4, 10.5), ylim=(0.58, 1.02),
           xlabel=r"$\beta_t$", ylabel=r"$A_1(\beta_t),\ A_2(\beta_t)$")
    ax.set_xticks([0.0, 2.5, 5.0, 7.5, 10.0])
    ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
    ax.xaxis.set_minor_locator(MultipleLocator(1.25))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_axisbelow(True)
    ax.grid(which="major", color="0.9", linewidth=1.0)
    ax.grid(which="minor", color="0.93", linewidth=0.6)
    ax.tick_params(which="minor", length=0)
    ax.legend(title="Acceptance rates", loc="center left",
              bbox_to_anchor=(1.03, 0.5), frameon=False,
              fontsize=11, labelspacing=1.5)

    output = Path(__file__).resolve().parent / "figures"
    output.mkdir(exist_ok=True)
    for extension in ("png", "svg"):
        path = output / f"Figure_003.{extension}"
        fig.savefig(path, dpi=300)
        print(f"Saved {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
