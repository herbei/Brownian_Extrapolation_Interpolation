 # Written by R Herbei (c) 2025
# to accompany Herbei & Somnath (2025)

  library(ggplot2)
  library(gridExtra)
  library(latex2exp)

  source("my_functions.R")

  alpha  <- seq(0.01, 10, length.out = 500)
  beta   <- pi^2 / (4 * alpha)
  alpha0 <- 1.091082

  s2 <- vapply(beta,  myseries2, numeric(1))
  s1 <- vapply(alpha, myseries1, numeric(1))

  II1 <- 4 * alpha * (1 - exp(-pi^2 / (4 * alpha))) / pi^2 + s2
  II2 <- 2 * pi^(-0.5) * alpha^(1.5) * s1

  df <- data.frame(alpha = alpha, z1 = II1, z2 = II2)

  p <- ggplot(df, aes(x = alpha)) +
    geom_line(aes(y = z1, colour = "z1"), linewidth = 0.5) +
    geom_line(aes(y = z2, colour = "z2"), linewidth = 0.5) +
    geom_point(aes(x = alpha0, y = 0), size = 2, colour = "red") +
    geom_vline(xintercept = alpha0, linetype = "dashed", colour = "grey") +
    annotate("text", x = alpha0, y = -0.1, label = "1.091082", vjust = 1.7, size = 3.5) +
    labs(x = TeX("$\\alpha_t$"), y = "", colour = "") +
    scale_colour_manual(
      labels = c(TeX("$\\frac{2I_1(\\alpha_t)}{a^2}$"), TeX("$\\frac{2I_2(\\alpha_t)}{a^2}$")),
      values = c("z1" = "red", "z2" = "blue")
    ) +
    coord_cartesian(ylim = c(-0.2, 2.5)) +
    theme_bw() +
    theme(
      legend.key.height = grid::unit(1.6, "cm"),
      legend.key.width  = grid::unit(3.5, "lines")
    )

  ggsave(
    filename = "Figures/Figure_002.png",
    plot     = p,
    width    = 6,
    height   = 3,
    units    = "in",
    dpi      = 300
  )
