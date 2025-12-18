library(ggplot2)
library(latex2exp)

alpha <- seq(0.001, 10, length.out = 500)
beta <- seq(0.001, 10, length.out = 500)
tol <- 1e-16

y_min <- min(alpha)
y_max <- max(alpha)
reverse_axis <- function(y) y_min + y_max - y

# Compute N for each curve; x-axis is log10(N)
alpha_log10N <- log10(sqrt((2 / alpha) * log(4 / (alpha * tol))))
beta_log10N <- log10(sqrt((1 / beta) * log(2 / tol)))

alpha_df <- data.frame(log10N = alpha_log10N, y = alpha, curve = "alpha")
beta_df <- data.frame(log10N = beta_log10N, y = reverse_axis(beta), curve = "beta")

curves <- rbind(alpha_df, beta_df)

combined_plot <- ggplot(curves, aes(x = log10N, y = y, color = curve, linetype = curve)) +
  geom_line(linewidth = 0.6) +
  labs(
    x = TeX("$\\log_{10}(N)$"),
    y = TeX("$\\alpha$"),
    color = NULL,
    linetype = NULL
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ reverse_axis(.), name = TeX("$\\beta$"))
  ) +
  scale_color_manual(
    values = c(alpha = "black", beta = "gray40"),
    labels = c(
      alpha = TeX("$\\log_{10}(N(\\alpha))$"),
      beta = TeX("$\\log_{10}(N(\\beta))$")
    )
  ) +
  scale_linetype_manual(
    values = c(alpha = "solid", beta = "dashed"),
    labels = c(
      alpha = TeX("$\\log_{10}(N(\\alpha))$"),
      beta = TeX("$\\log_{10}(N(\\beta))$")
    )
  ) +
  theme_bw()

ggsave(
  filename = "Figures/Figure_009.png",
  plot = combined_plot,
  width = 6,
  height = 2.5,
  units = "in",
  dpi = 300
)
