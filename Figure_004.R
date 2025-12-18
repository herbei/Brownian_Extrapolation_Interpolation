rm(list = ls())

# Work from the script directory so relative paths resolve when run from anywhere
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(script_dir)

library(ggplot2)
library(latex2exp)

source('my_functions.R')

simulate_draws <- function(t, a, B) {
  draws <- vapply(seq_len(B), function(i) extrapolate_bm_II(t, a), numeric(1))
  data.frame(y = draws, t = t)
}

theoretical_density <- function(t, a, grid_n = 200) {
  z_vals <- seq(0, a, length.out = grid_n)
  dz <- z_vals[2] - z_vals[1]
  d_vals <- h_series(z_vals, t, a)
  norm_constant <- dz * sum(d_vals)
  data.frame(z = z_vals, Z = d_vals / norm_constant, t = t)
}

a <- 2
B <- 10000
t_vals <- c(0.2, 1)

sample_data <- do.call(rbind, lapply(t_vals, simulate_draws, a = a, B = B))
density_data <- do.call(rbind, lapply(t_vals, theoretical_density, a = a))

facet_labeller <- labeller(t = function(x) paste0("t = ", x))

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
  labs(x = "y", y = TeX("$h(y;t,a)$")) +
  coord_cartesian(xlim = c(0, a)) +
  facet_wrap(~t, nrow = 1, labeller = facet_labeller) +
  theme_bw()

print(p)

ggsave(
  filename = "Figures/Figure_004.png",
  plot = p,
  width = 8,
  height = 3,
  units = "in",
  dpi = 300
)
