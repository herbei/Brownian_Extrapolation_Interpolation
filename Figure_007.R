rm(list = ls())

# Work from the script directory so relative paths resolve when run from anywhere
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(script_dir)

library(patchwork)
library(ggplot2)

# Compatibility helper for older ggplot2 (before after_stat was available)
if (!exists("after_stat", mode = "function")) {
  after_stat <- function(x) x
}

p_a <- readRDS("Figures/p_007a.rds")
p_b <- readRDS("Figures/p_007b.rds")

p <- p_a + p_b + plot_layout(nrow = 1, widths = c(1, 1.1))

print(p)

ggsave(
  filename = "Figures/Figure_007.png",
  plot = p,
  width = 14,
  height = 3,
  units = "in",
  dpi = 300
)
