# Rcode README (Figures + Monte Carlo)

This folder contains:
- `my_functions.R`: core functions used by all scripts (do not edit unless you intend to change the methods).
- `Figure_*.R`: scripts that generate the paper figures into `Figures/`.
- `mc_interpolate_BB_*.R`: Monte Carlo scripts that benchmark the rejection samplers used by `interpolate_BB_*`.

## Requirements

- R (recent version recommended)
- Packages used across scripts:
  - `ggplot2`, `latex2exp`
  - `patchwork` (used by `Figure_005.R`, `Figure_006.R`, `Figure_007.R`, `Figure_008.R`)
  - `gridExtra` (used by `Figure_002.R`)

Install packages in R if needed:
```r
install.packages(c("ggplot2", "latex2exp", "patchwork", "gridExtra"))
```

## Figure scripts (`Figure_*.R`)

Run from this directory (so relative paths like `Figures/...` resolve):
```sh
Rscript Figure_001.R
```

To (re)generate all figures:
```sh
for f in Figure_*.R; do Rscript "$f"; done
```

**Outputs**
- Most scripts write `.png` files to `Figures/` via `ggsave()`.
- Some scripts also save plot objects as `.rds` in `Figures/` (e.g., `Figures/p_005a.rds`).

## Monte Carlo scripts (`mc_interpolate_BB_*.R`)

These scripts report, for each scenario:
1. **Acceptance rate** of the rejection sampler inside `interpolate_BB_*`, estimated as:
   - `acceptance_rate = B_accept / total_proposals`
   - where `total_proposals` is counted by temporarily wrapping the proposal generator (`extrapolate_bm_I` or `extrapolate_bm_II`).
2. **Elapsed wall-clock time** to generate `B_time = 10000` accepted draws, using `system.time(...)[["elapsed"]]`.

Run individually:
```sh
Rscript mc_interpolate_BB_I.R
Rscript mc_interpolate_BB_II.R
Rscript mc_interpolate_BB_III.R
Rscript mc_interpolate_BB_IV.R
```

**What varies in each script**
- `mc_interpolate_BB_I.R`: fixes `a=2`, `T=3`, uses `t ∈ {0.2, 1.5, 2.8}`, varies `x` and `z` over grids in `(0,a)`.
- `mc_interpolate_BB_II.R`: fixes `a=2`, `T=3`, uses `t ∈ {0.2, 1.5, 2.8}`, varies `z` over a grid in `(0,a)`.
- `mc_interpolate_BB_III.R`: fixes `a=2`, `T=3`, uses `t ∈ {0.2, 1.5, 2.8}`, varies `x` over a grid in `(0,a)`.
- `mc_interpolate_BB_IV.R`: fixes `a=2`, `T=3`, uses `t ∈ {0.2, 1.5, 2.8}` (no endpoint grid in this case).

**Outputs**
Each script writes two CSV files with a timestamped `run_id`:
- `mc_interpolate_BB_*.csv`: per-scenario results (includes acceptance rate and elapsed seconds for 10,000 draws).
- `mc_interpolate_BB_*_summary_by_t_<run_id>.csv`: a simple summary by `t`.
  - For `mc_interpolate_BB_I.R`, the summary is averaged over the `(x,z)` grid:
    - mean acceptance rate over `(x,z)`
    - median elapsed time (over `(x,z)`) to generate 10,000 draws

## Customizing the experiments

In each `mc_interpolate_BB_*.R` file you can edit:
- grid sizes (`x_grid_n`, `z_grid_n`) and the boundary buffer `eps`
- `B_accept` (acceptance-rate accuracy)
- `B_time` (timed number of accepted draws; currently 10,000)
- `timing_reps` (set `>1` to repeat timings and take the median across repeats)

