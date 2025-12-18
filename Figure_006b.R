rm(list = ls())

# Work from the script directory so relative paths resolve when run from anywhere
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(script_dir)

set.seed(4)

library(ggplot2)
library(latex2exp)

# Compatibility helper for older ggplot2 (before after_stat was available)
if (!exists("after_stat", mode = "function")) {
  after_stat <- function(x) x
}

source("my_functions.R")

simulate_draws <- function(t, T, z, a, B) {
  draws <- vapply(seq_len(B), function(i) interpolate_BB_II(t, T, z, a), numeric(1))
  data.frame(y = draws, t = t)
}

theoretical_overlay <- function(t, T, z, a, w, grid_n = 200) {
  z_vals <- seq(0, a, length.out = grid_n)
  dz <- z_vals[2] - z_vals[1]
  d_vals <- h_series(z_vals, t, a) * g_series(z_vals, T - t, z, a)
  norm_constant <- dz * sum(d_vals)
  data.frame(x = t + w * (d_vals / norm_constant), y = z_vals, t = t)
}

rect_df <- function(z, t, bins, w) {
  h <- hist(z, breaks = bins, plot = FALSE)
  data.frame(
    xmin = t,
    xmax = t + w * h$density,
    ymin = h$breaks[-length(h$breaks)],
    ymax = h$breaks[-1],
    t = t
  )
}

kde_df <- function(z, t, w) {
  d <- density(z)
  data.frame(
    x = t + w * d$y,
    y = d$x,
    t = t
  )
}

pick_points <- function(path_df, times) {
  idx <- vapply(times, function(tv) which.min(abs(path_df$t_vals - tv)), integer(1))
  data.frame(x = path_df$t_vals[idx], y = path_df$BB_vals[idx])
}

T_total <- 3
a <- 2
z <- 1.8
# Offsets for left/right annotations to keep text off the axes when paneled
left_label_x <- -0.35
right_label_x <- T_total + 0.3
t_vals <- c(0.2, 2.8)
B <- 10000
w <- 0.5
bins <- 30

sample_data <- do.call(rbind, lapply(t_vals, simulate_draws, T = T_total, z = z, a = a, B = B))

rects <- do.call(
  rbind,
  lapply(split(sample_data, sample_data$t), function(df) rect_df(df$y, unique(df$t), bins, w))
)

kdes <- do.call(
  rbind,
  lapply(split(sample_data, sample_data$t), function(df) kde_df(df$y, unique(df$t), w))
)

overlay_data <- do.call(
  rbind,
  lapply(t_vals, theoretical_overlay, T = T_total, z = z, a = a, w = w)
)

paths <- replicate(5, sim_bb_path(0, z, T_total, a, N = 300), simplify = FALSE)
paths_df <- do.call(rbind, Map(function(df, id) transform(df, path = id), paths, seq_along(paths)))
points_df <- do.call(
  rbind,
  Map(function(df, id) transform(pick_points(df, t_vals), path = id), paths, seq_along(paths))
)

boundary_segments <- data.frame(
  x = c(0, 0, T_total, 0),
  y = c(0, 0, 0, a),
  xend = c(0, T_total, T_total, T_total),
  yend = c(a, 0, a, a)
)

endpoint_data <- data.frame(
  x = c(T_total, 0),
  y = c(z, 0)
)

my_alt_colors <- c("teal" = "#1b9e77", "goldenrod" = "#d95f02", "skyblue" = "#7570b3",
                   "crimson" = "#e7298a", "olive" = "#66a61e")
path_colors <- setNames(unname(my_alt_colors), as.character(seq_along(my_alt_colors)))
paths_df$path <- factor(paths_df$path)
points_df$path <- factor(points_df$path)

p <- ggplot() +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = t),
    fill = "grey75",
    color = "grey35",
    alpha = 0.3
  ) +
  geom_path(data = kdes, aes(x = x, y = y, group = t), linewidth = 0.5, color = "blue") +
  geom_path(data = overlay_data, aes(x = x, y = y, group = t), color = "red") +
  geom_line(data = paths_df, aes(x = t_vals, y = BB_vals, color = factor(path)), linewidth = 0.3, alpha = 0.5) +
  geom_point(data = points_df, aes(x = x, y = y, color = factor(path)), size = 2, show.legend = FALSE) +
  geom_point(data = endpoint_data, aes(x = x, y = y), size = 2, color = "blue", inherit.aes = FALSE) +
  geom_segment(data = boundary_segments, aes(x = x, y = y, xend = xend, yend = yend), color = "darkgreen") +
  geom_vline(xintercept = t_vals, linewidth = 0.3, linetype = 2) +
  labs(x = "t", y = "W(t)") +
  coord_cartesian(xlim = c(left_label_x, T_total + 0.6), clip = "off") +
  annotate("text", x = right_label_x, y = 0.1, label = TeX("$T=3$"), color = "red", parse = TRUE, size = 5) +
  annotate("text", x = left_label_x, y = a + 0.1, label = TeX("$a=2$"), color = "red", parse = TRUE, size = 5) +
  annotate("text", x = right_label_x, y = z + 0.1, label = TeX("$z=1.8$"), color = "red", parse = TRUE, size = 5) +
  annotate("text", x = left_label_x, y = 0.1, label = TeX("$x=0$"), color = "red", parse = TRUE, size = 5) +
  annotate("text", x = t_vals[1], y = -0.1, label = TeX("$t_1 = 0.2$"), color = "red", parse = TRUE, size = 5) +
  annotate("text", x = t_vals[2], y = -0.1, label = TeX("$t_2 = 2.8$"), color = "red", parse = TRUE, size = 5) +
  scale_color_manual(values = path_colors, guide = "none") +
  theme_bw()

print(p)

ggsave(
  filename = "Figures/Figure_006b.png",
  plot = p,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

saveRDS(p, "Figures/p_006b.rds")
