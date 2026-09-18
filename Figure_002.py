"""Reproduce Figure 2 and calculate alpha_0 and beta_0.

Run: python3 Figure_002.py
Dependencies: numpy, scipy, matplotlib.

The curves are J1 = 2*I1/a**2 and J2 = 2*I2/a**2 from Section 3.1.
They depend only on alpha, with beta = pi**2/(4*alpha).
The 200-term sums are sufficient for the plotted range 0.01 <= alpha <= 10;
this script is a numerical figure calculation, not a sampling routine.
"""

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from scipy.optimize import brentq


def scaled_integrals(alpha, terms=200):
    """Return the scaled envelope integrals (2*I1/a**2, 2*I2/a**2)."""
    beta = math.pi**2 / (4.0 * alpha)
    s1 = math.fsum(n**2 * math.exp(-alpha * n**2)
                   for n in range(1, terms + 1))
    s2 = math.fsum(math.exp(-beta * n**2)
                   for n in range(1, terms + 1))
    j1 = -math.expm1(-beta) / beta + s2
    j2 = 2.0 / math.sqrt(math.pi) * alpha**1.5 * s1
    return j1, j2


def calculate_cutoffs(terms=200):
    """Solve J1(alpha_0) = J2(alpha_0), then transform to beta_0."""
    def difference(alpha):
        j1, j2 = scaled_integrals(alpha, terms)
        return j1 - j2

    alpha0 = brentq(difference, 1.0, 1.2, xtol=1e-14)
    beta0 = math.pi**2 / (4.0 * alpha0)
    return alpha0, beta0


def main():
    alpha0, beta0 = calculate_cutoffs()
    j1_at_root, j2_at_root = scaled_integrals(alpha0)
    print(f"alpha_0 = {alpha0:.12f}")
    print(f"beta_0  = {beta0:.12f}")
    print(f"Common scaled integral = {j1_at_root:.12f}")
    print(f"J1(alpha_0) - J2(alpha_0) = {j1_at_root - j2_at_root:.3e}")

    alpha = np.linspace(0.01, 10.0, 500)
    curves = np.array([scaled_integrals(value) for value in alpha])

    fig, ax = plt.subplots(figsize=(6, 3))
    fig.subplots_adjust(left=0.08, right=0.74, bottom=0.20, top=0.96)
    ax.plot(alpha, curves[:, 0], color="red", linewidth=1.3,
            label=r"$\frac{2I_1(\alpha_t)}{a^2}$")
    ax.plot(alpha, curves[:, 1], color="blue", linewidth=1.3,
            label=r"$\frac{2I_2(\alpha_t)}{a^2}$")
    ax.axvline(alpha0, color="0.7", linestyle="--", linewidth=1.2)
    ax.plot(alpha0, 0.0, "o", color="red", markersize=5)
    ax.text(alpha0, -0.17, f"{alpha0:.6f}", ha="center", va="top")
    ax.set(xlim=(-0.45, 10.5), ylim=(-0.33, 2.63), xlabel=r"$\alpha_t$")
    ax.set_xticks([0.0, 2.5, 5.0, 7.5, 10.0])
    ax.set_yticks([0.0, 1.0, 2.0])
    ax.xaxis.set_minor_locator(MultipleLocator(1.25))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_axisbelow(True)
    ax.grid(which="major", color="0.9", linewidth=1.0)
    ax.grid(which="minor", color="0.93", linewidth=0.6)
    ax.tick_params(which="minor", length=0)
    ax.legend(loc="center left", bbox_to_anchor=(1.03, 0.5),
              frameon=False, fontsize=12, labelspacing=2.5, handlelength=2.5)

    output = Path(__file__).resolve().parent / "figures"
    output.mkdir(exist_ok=True)
    for extension in ("png", "svg"):
        path = output / f"Figure_002.{extension}"
        fig.savefig(path, dpi=300)
        print(f"Saved {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
