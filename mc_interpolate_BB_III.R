rm(list = ls())

# Work from the script directory so relative paths resolve when run from anywhere
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(script_dir)

set.seed(1)

source("my_functions.R")

# ----------------------------
# Experiment configuration
# ----------------------------

# Fixed parameters (per your spec)
a_fixed <- 2
T_fixed <- 3
t_vals <- c(0.2, 1.5, 2.8)

# Grid over x in (0, a)
eps <- 0.05
x_grid_n <- 5L
x_grid <- seq(eps, a_fixed - eps, length.out = x_grid_n)

# Number of accepted draws for acceptance-rate estimation
# (does not need to match B_time; increase if you want tighter MC error)
B_accept <- 2000L

# Number of accepted draws for timing (requested: 10000)
B_time <- 10000L

# Repetitions for timing to reduce noise (report median)
timing_reps <- 1L

scenarios <- expand.grid(
  t = t_vals,
  x = x_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
scenarios$T <- T_fixed
scenarios$a <- a_fixed
scenarios$scenario_id <- seq_len(nrow(scenarios))
scenarios <- scenarios[, c("scenario_id", "t", "T", "x", "a")]

# ----------------------------
# Helpers
# ----------------------------

count_proposals_per_draw <- function(t, T, x, a, B) {
  proposals_env <- new.env(parent = emptyenv())
  proposals_env$n <- 0L

  base_extrapolate <- get("extrapolate_bm_II", envir = .GlobalEnv)
  assign(
    "extrapolate_bm_II",
    function(...) {
      proposals_env$n <- proposals_env$n + 1L
      base_extrapolate(...)
    },
    envir = .GlobalEnv
  )
  on.exit(assign("extrapolate_bm_II", base_extrapolate, envir = .GlobalEnv), add = TRUE)

  proposals_per_draw <- integer(B)
  for (i in seq_len(B)) {
    before <- proposals_env$n
    interpolate_BB_III(t, T, x, a)
    proposals_per_draw[i] <- proposals_env$n - before
  }

  proposals_per_draw
}

time_draws_elapsed <- function(t, T, x, a, B, reps = 1L, scenario_id = 1L) {
  elapsed <- numeric(reps)
  for (r in seq_len(reps)) {
    set.seed(200000L + as.integer(scenario_id) * 100L + r)
    invisible(gc())
    out <- numeric(B)
    elapsed[r] <- system.time({
      for (i in seq_len(B)) {
        out[i] <- interpolate_BB_III(t, T, x, a)
      }
    })[["elapsed"]]
  }
  elapsed
}

# ----------------------------
# Run experiment
# ----------------------------

results <- vector("list", nrow(scenarios))

for (k in seq_len(nrow(scenarios))) {
  sc <- scenarios[k, ]
  set.seed(100000L + sc$scenario_id)
  cat(
    sprintf(
      "Scenario %d/%d: t=%.3f, T=%.3f, x=%.3f, a=%.3f\n",
      k, nrow(scenarios), sc$t, sc$T, sc$x, sc$a
    )
  )

  proposals_per_draw <- count_proposals_per_draw(sc$t, sc$T, sc$x, sc$a, B_accept)
  total_proposals <- sum(proposals_per_draw)
  acceptance_rate <- B_accept / total_proposals

  elapsed <- time_draws_elapsed(
    sc$t, sc$T, sc$x, sc$a,
    B_time,
    reps = timing_reps,
    scenario_id = sc$scenario_id
  )

  results[[k]] <- data.frame(
    scenario_id = sc$scenario_id,
    t = sc$t,
    T = sc$T,
    x = sc$x,
    a = sc$a,
    B_accept = B_accept,
    total_proposals = total_proposals,
    acceptance_rate = acceptance_rate,
    mean_proposals_per_draw = mean(proposals_per_draw),
    median_proposals_per_draw = stats::median(proposals_per_draw),
    p90_proposals_per_draw = stats::quantile(proposals_per_draw, 0.9, names = FALSE, type = 7),
    B_time = B_time,
    timing_reps = timing_reps,
    elapsed_median_sec = stats::median(elapsed),
    elapsed_mean_sec = mean(elapsed),
    elapsed_sd_sec = stats::sd(elapsed)
  )
}

results_df <- do.call(rbind, results)
print(results_df)

summary_by_t <- do.call(
  rbind,
  lapply(sort(unique(results_df$t)), function(tv) {
    df <- results_df[results_df$t == tv, , drop = FALSE]
    data.frame(
      t = tv,
      acceptance_rate_mean_x = mean(df$acceptance_rate),
      elapsed_median_sec_median_x = stats::median(df$elapsed_median_sec),
      stringsAsFactors = FALSE
    )
  })
)

cat("\nSummary by t (averaged over x grid):\n")
print(summary_by_t, row.names = FALSE)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_csv <- sprintf("mc_interpolate_BB_III_results_%s.csv", run_id)
write.csv(results_df, out_csv, row.names = FALSE)
cat(sprintf("Wrote results to %s\n", out_csv))

summary_csv <- sprintf("mc_interpolate_BB_III_summary_by_t_%s.csv", run_id)
write.csv(summary_by_t, summary_csv, row.names = FALSE)
cat(sprintf("Wrote summary to %s\n", summary_csv))

