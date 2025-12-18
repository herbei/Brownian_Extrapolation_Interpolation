rm(list = ls())

# Work from the script directory so relative paths resolve when run from anywhere
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(script_dir)

set.seed(1)

library(ggplot2)
library(latex2exp)

# Compatibility helper for older ggplot2 (before after_stat was available)
if (!exists("after_stat", mode = "function")) {
  after_stat <- function(x) x
}

source("my_functions.R")

simulate_draws <- function(t, T, x, z, a, B) {
  draws <- vapply(seq_len(B), function(i) interpolate_BB_I(t, T, x, z, a), numeric(1))
  data.frame(y = draws, t = t)
}

theoretical_density <- function(t, T, x, z, a, grid_n = 200) {
  z_vals <- seq(0, a, length.out = grid_n)
  dz <- z_vals[2] - z_vals[1]
  d_vals <- g_series(z_vals, t, x, a) * g_series(z_vals, T - t, z, a)
  norm_constant <- dz * sum(d_vals)
  data.frame(z = z_vals, Z = d_vals / norm_constant, t = t)
}

T <- 3
x <- 0.1
z <- 1.8
a <- 2
B <- 10000
t_vals <- c(0.2, 2.8)

sample_data <- do.call(rbind, lapply(t_vals, simulate_draws, T = T, x = x, z = z, a = a, B = B))
density_data <- do.call(rbind, lapply(t_vals, theoretical_density, T = T, x = x, z = z, a = a))

facet_labeller <- labeller(t = function(val) paste0("t = ", val))

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
  labs(x = "y", y = TeX("$k_{(1)}(y)$")) +
  coord_cartesian(xlim = c(0, a)) +
  facet_wrap(~t, nrow = 1, labeller = facet_labeller) +
  theme_bw()

print(p)

ggsave(
  filename = "Figures/Figure_005a.png",
  plot = p,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

saveRDS(p, "Figures/p_005a.rds")
