# Interval-Constrained Brownian Paths

Python implementation accompanying **Interval-Constrained Brownian Paths:
Exact Interpolation and Extrapolation**, by Radu Herbei and Kumar Somnath.

The code samples Brownian motion endpoints conditional on survival in a fixed
interval, and Brownian bridge values conditional on remaining in that interval.
It also reproduces Figures 1–9 and the acceptance-rate and timing experiment
in Section 3.3. The scripts use ordinary Python functions and loops.

## Files

| File or directory | Purpose |
| --- | --- |
| `ICB_gh_samplers.py` | Extrapolation samplers, image/spectral comparison bounds, and one-step bridge interpolation (Algorithms 1–6). |
| `ICB_samplers.py` | Sequential bridge skeletons on a general interval `[l,u]`, plus numerical kernel and unconstrained-bridge helpers. |
| `Figure_001.py` through `Figure_009.py` | One standalone script for each manuscript figure. |
| `timing_experiment.py` | Section 3.3 acceptance estimates and timings, with pilot and full modes. |
| `tests/` | Regression checks for endpoints, distributions, image bounds, numerical behavior, time reversal, and the timing procedure. |
| `requirements.txt` | Third-party packages needed for the figures. |
| `figures/` | Reference PNG/SVG figures; `figures/data/` contains the generated CSV data. |
| `timing_results/full_20260917_121103_136179/` | The full experiment reported in the manuscript, including measured-source snapshots. |
| `MANIFEST.csv`, `SHA256SUMS.txt` | File inventory and checksums for the prepared release. |

The release is a collection of adjacent scripts; no package installation is
needed. Keep the two sampler modules beside the figure and timing scripts.

## Requirements

Use **Python 3.9 or later**. The samplers, timing experiment, and tests use only
the Python standard library. Figure generation additionally needs NumPy,
SciPy, and Matplotlib:

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\activate` and use
`python` in place of `python3` if appropriate.

The reference results were generated with Python 3.9.6, NumPy 2.0.2,
SciPy 1.13.1, and Matplotlib 3.9.4 on macOS/arm64. These are the tested
versions; `requirements.txt` allows the packages available for your Python
version. Plot rendering and numerical results can differ across environments.

## Sampling examples

Run these examples from the directory containing the scripts:

```python
import random
from ICB_gh_samplers import ICB_sample_g, ICB_sample_h, ICB_interpolate_step
from ICB_samplers import ICB

random.seed(12345)

# Interior-start extrapolation on [0, 2].
y_g = ICB_sample_g(t=0.2, x=0.5, a=2.0)

# Lower-boundary-start extrapolation on [0, 2].
y_h = ICB_sample_h(t=0.2, a=2.0)

# One bridge value at time 1, given endpoints 0.1 and 1.8 at times 0 and 3.
y_bridge = ICB_interpolate_step(t=1.0, T=3.0, x=0.1, z=1.8, a=2.0)

# Two independent bridge skeletons, each with three interior observations.
rows = ICB(
    Delta=3.0, x=0.1, z=1.8, l=0.0, u=2.0,
    calT=[0.5, 1.0, 2.5], n=2,
)
assert len(rows) == 10  # Three interior times plus both endpoints per draw.
print(rows[0])
```

`ICB` returns a list of dictionaries with keys `draw`, `time`, `value`, and
`attempts`, ordered by draw and then time. Draw numbers start at one. The
two endpoints are included by default; set `include_endpoints=False` to omit
them. The requested interior times must be distinct and strictly between
zero and `Delta`; a sorted copy is used internally. The `attempts` field is
always one and does not count rejection proposals.

The common one-step interface `ICB_interpolate_step` dispatches to the
appropriate endpoint case, using reflection and time reversal as needed:

| Manuscript algorithm | Direct function | Endpoint case on `[0,a]` |
| --- | --- | --- |
| 3 | `ICB_interpolate_BB_I(t, T, x, z, a)` | Both endpoints interior. |
| 4 | `ICB_interpolate_BB_II(t, T, z, a)` | Lower boundary to interior. |
| 5 | `ICB_interpolate_BB_III(t, T, x, a)` | Interior to upper boundary. |
| 6 | `ICB_interpolate_BB_IV(t, T, a)` | Lower boundary to upper boundary. |

## Implementation choices

The implementations of Algorithms 3–6 enable the time-reversal option
described in Section 3.3, selecting the proposal orientation associated
with the shorter time interval. The figure simulations and timing results
use this option. The manuscript displays the basic algorithms for simplicity.

For rejection comparisons involving `g` and `h`, the code uses the image
bounds in Propositions 6–7 when `alpha = pi**2*t/(2*a**2) < pi`, and the
spectral bounds in Propositions 4–5 when `alpha >= pi`. Each factor uses its
own time argument. The image interval is `S_N - error_N` to `S_N + error_N`
for `g`, and `S_N - error_N` to `S_N` for `h`. An unresolved comparison
increases the truncation index by one.

This representation switch is separate from the boundary-start proposal
cutoffs `alpha_0 = 1.0861647656269326` and
`beta_0 = pi**2/(4*alpha_0) = 2.271663727600443`, calculated by
`Figure_002.py`. The interior-start proposal uses the regimes
`a/sqrt(t) < 2` and `a/sqrt(t) >= 2`.

The sampling algorithms are mathematically exact series-method rejection
procedures, implemented using ordinary floating-point arithmetic. Numerical
plotting helpers use finite sums and are not exact density evaluators.
In particular, `ICB_g_series` and `ICB_h_series` evaluate finite spectral
sums; the image-series helpers are separate. Where supported, `scaled=True`
multiplies the kernel or bound by `exp(alpha)` to avoid underflow. It does
not normalize a density.

Same-boundary bridge interpolation (`0` to `0`, or `a` to `a`) is not
implemented and raises `NotImplementedError` when interior observations
are requested. Only exact equality identifies a boundary endpoint.
Extreme floating-point inputs can exceed the supported numerical range;
several underflow or loss-of-interior-position cases explicitly raise
`FloatingPointError`. In `ICB`, `terms` and `rel_tol` are compatibility
parameters and do not control the exact sampler; `max_attempts` must remain
positive infinity.

## Reproduce Figures 1–9

Run each script separately, or run all nine in order:

```sh
python3 Figure_001.py
python3 Figure_002.py
python3 Figure_003.py
python3 Figure_004.py
python3 Figure_005.py
python3 Figure_006.py
python3 Figure_007.py
python3 Figure_008.py
python3 Figure_009.py
```

Each script writes PNG and SVG files under `figures/`, relative to the
script location. Rerunning a script replaces its corresponding outputs.
Figures 1 and 4–8 also save samples and plotted density values as CSV;
Figures 5–8 additionally save path values. Figure 9 saves its bound curves
and integer ceilings. Figures 2–3 are deterministic and do not write CSVs.
The release includes all 17 generated CSV files.

| Figure | Default setup | Sampling seed |
| --- | --- | --- |
| 1 | `x=0.5`; left `(a,t)=(1.8,1)`; right `(2.5,0.2)`. | 1001 |
| 2 | Envelope integrals and the proposal cutoffs `alpha_0`, `beta_0`. | Deterministic |
| 3 | Theoretical acceptance curves for Algorithms 1 and 2. | Deterministic |
| 4 | `a=2`, with `t=0.2` and `t=1`. | 1004 |
| 5 | `a=2`, `T=3`, `x=0.1`, `z=1.8`. | 2026091605 |
| 6 | `a=2`, `T=3`, `x=0`, `z=1.8`. | 2026091606 |
| 7 | `a=2`, `T=3`, `x=0.5`, `z=2`. | 7001 |
| 8 | `a=2`, `T=3`, `x=0`, `z=2`. | 8001 |
| 9 | Continuous truncation expressions with `tol=1e-16`. | Deterministic |

Figures 1 and 4 use 10,000 draws per panel. Figures 5–8 use 10,000 joint
pairs at times 0.2 and 2.8, and five separately sampled bridge skeletons
on 301 times including both endpoints. The path seed is the sampling seed
plus 1000. Lines joining skeleton values are for display. Parameters and
seeds are near the top of each script or in its `main()` function.

## Reproduce the Section 3.3 experiment

For a pilot, with 100 accepted draws for acceptance estimation and a
separate batch of 500 accepted draws for timing per setting:

```sh
python3 timing_experiment.py --pilot
```

For the full experiment, with 2,000 accepted draws for acceptance estimation
and a separate batch of 10,000 accepted draws for timing per setting:

```sh
python3 timing_experiment.py
```

Both modes use `a=2`, `T=3`, the time grid
`0.1, 1, 1.5, 2, 2.5, 2.9`, and the interior endpoint grid
`0.05, 0.525, 1, 1.475, 1.95`. There are 216 settings: 150 for Algorithm 3,
30 each for Algorithms 4 and 5, and six for Algorithm 6. One timing batch
is measured per setting. The full run generates 432,000 acceptance draws
and 2,160,000 timed draws, excluding rejected proposals. It took about
30 minutes in the recorded environment.

Acceptance is the number of accepted draws divided by the number of outer
bridge proposals; rejections internal to `g` or `h` sampling are excluded.
The table averages the per-setting acceptance ratios over endpoint settings.
Timing uses `time.perf_counter()` around the sampling-and-storage loop,
with proposal counting disabled. Seeding, output allocation, validation,
and file writing are outside the timed loop. The table takes the median
batch time across endpoint settings. Times are seconds per 10,000 draws
in full mode and per 500 draws in pilot mode.

Runs are serial and create fresh UTC-stamped folders under `timing_results/`
containing `scenarios.csv`, `timing_batches.csv`, `summary.csv`, `summary.md`,
and `metadata.json`. They do not edit the manuscript. Acceptance seeds are
`100000 + scenario_id`; timing seeds are
`200000 + 100*scenario_id + repetition`, starting with repetition one.
Scenario IDs restart for each algorithm, with time varying fastest, then
`x`, then `z`. Timing results depend on hardware and system load.

The recorded full results supporting the manuscript are in
[`timing_results/full_20260917_121103_136179/`](timing_results/full_20260917_121103_136179/summary.md).
Its `source_snapshot/` retains the exact two source files measured in that
run; their hashes match `metadata.json`. The current source differs only
in documentation, and the bridge-wrapper import was made usable directly
from this directory. The release copy of the metadata replaces the original
absolute executable path with `python3` and records this change explicitly.
The original benchmark measurements and source hashes are preserved.

## Run the tests

```sh
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

The suite contains 40 tests and requires no third-party packages. It checks
finite-precision edge cases, endpoint handling, the comparison bounds,
time-reversal behavior, representative distributional properties, and the
benchmark's counting/summary logic. These checks support reproducibility
and regression detection; the mathematical arguments are in the manuscript.

## Reference

Radu Herbei and Kumar Somnath (2026). *Interval-Constrained Brownian Paths:
Exact Interpolation and Extrapolation*. Manuscript.

The interior-start proposal draws on Luc Devroye (2010), *On Exact Simulation
Algorithms for Some Distributions Related to Brownian Motion and Brownian
Meanders*, DOI: [10.1007/978-3-7908-2598-5_1](https://doi.org/10.1007/978-3-7908-2598-5_1).
