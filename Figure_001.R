rm(list = ls())

# Work from the script directory so relative paths resolve when run from anywhere
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(script_dir)

set.seed(1)

library(ggplot2)
library(latex2exp)

source("my_functions.R")

simulate_draws <- function(a, t, x, B) {
  draws <- vapply(seq_len(B), function(i) extrapolate_bm_I(t, x, a), numeric(1))
  data.frame(y = draws, a = a)
}

theoretical_density <- function(a, t, x, grid_n = 200) {
  z_vals <- seq(0, a, length.out = grid_n)
  dz <- z_vals[2] - z_vals[1]
  d_vals <- g_series(z_vals, t, x, a)
  norm_constant <- dz * sum(d_vals)
  data.frame(z = z_vals, Z = d_vals / norm_constant, a = a)
}

t <- 0.2
x <- 0.5
B <- 10000
a_vals <- c(1.8, 2.5)

sample_data <- do.call(rbind, lapply(a_vals, simulate_draws, t = t, x = x, B = B))
density_data <- do.call(rbind, lapply(a_vals, theoretical_density, t = t, x = x))

facet_labeller <- labeller(a = function(val) paste0("a = ", val))

p <- ggplot(sample_data, aes(x = y)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "steelblue",
    colour = "black",
    linewidth = 0.1,
    alpha = 0.5
  ) +
  geom_density(aes(y = after_stat(density)), linewidth = 0.5, color = "blue") +
  geom_line(
    data = density_data,
    aes(x = z, y = Z),
    linewidth = 0.5,
    color = "red",
    inherit.aes = FALSE
  ) +
  labs(x = "y", y = TeX("$g(y;t,x,a)$")) +
  facet_wrap(~a, nrow = 1, labeller = facet_labeller, scales = "free") +
  theme_bw()

print(p)

ggsave(
  filename = "Figures/Figure_001.png",
  plot = p,
  width = 8,
  height = 3,
  units = "in",
  dpi = 300
)
