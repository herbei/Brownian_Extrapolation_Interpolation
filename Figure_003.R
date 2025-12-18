rm(list = ls())

# Work from the script directory so relative paths resolve when run from anywhere
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(script_dir)
library(ggplot2)
library(latex2exp)

source('my_functions.R')

N <- 300
alpha0 <- 1.091082
beta0 <- pi^2 / (4 * alpha0)

series_numerator <- function(beta_vals) {
  vapply(beta_vals, myseries3, numeric(1))
}

denominator_A1 <- function(alpha_vals, beta_vals) {
  coeff <- (pi^(5 / 2) / 4) * beta_vals^(-1 / 2)
  coeff * vapply(alpha_vals, myseries1, numeric(1))
}

denominator_A2 <- function(alpha_vals, beta_vals) {
  1 - exp(-beta_vals) + beta_vals * vapply(beta_vals, myseries2, numeric(1))
}

calc_ACC <- function(beta_vals, denom_fn) {
  alpha_vals <- (pi^2) / (4 * beta_vals)
  numerators <- series_numerator(beta_vals)
  numerators / denom_fn(alpha_vals, beta_vals)
}

beta_vals <- seq(0.1, beta0, length.out = N)
beta_vals1 <- seq(beta0, 10, length.out = N)

ACC <- calc_ACC(beta_vals, denominator_A1)
ACC1 <- calc_ACC(beta_vals1, denominator_A2)

plot_data <- rbind(
  data.frame(beta = beta_vals, A = ACC, series = "A1"),
  data.frame(beta = beta_vals1, A = ACC1, series = "A2")
)

p3 <- ggplot(plot_data, aes(x = beta, y = A, color = series)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(
    values = c(A1 = "red", A2 = "blue"),
    labels = c(A1 = TeX("$A_2(\\beta_t)$"), A2 = TeX("$A_1(\\beta_t)$")),
    name = TeX("Acceptance rates")
  ) +
  geom_vline(xintercept = 2.261426, linewidth = 0.3, linetype = 2, color = "grey") +
  xlab(TeX("$\\beta_t$")) +
  ylab(TeX("$A_1(\\beta_t),\\ A_2(\\beta_t)$")) +
  geom_point(aes(x = beta0, y = 0.7), color = "red", inherit.aes = FALSE) +
  annotate("text", x = beta0, y = 0.6, label = TeX("$\\beta_0=2.26142$"), color = "black", parse = TRUE, size = 4) +
  annotate("text", x = 1, y = 0.7, label = TeX("$\\beta_t<\\beta_0$"), color = "red", parse = TRUE, size = 5) +
  annotate("text", x = 5, y = 0.7, label = TeX("$\\beta_t>\\beta_0$"), color = "blue", parse = TRUE, size = 5) +
  theme_bw()
print(p3)

ggsave(
  filename = "Figures/Figure_003.png",
  plot = p3,
  width = 5,
  height = 2,
  units = "in",       # or "cm", "mm"
  dpi = 300
)
