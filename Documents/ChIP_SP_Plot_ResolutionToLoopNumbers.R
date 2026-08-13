suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# outputs
out_png <- "HiC_loop_counts_5_10_25kb_by_cellline_smooth.png"
out_pdf <- "HiC_loop_counts_5_10_25kb_by_cellline_smooth.pdf"

summary_df <- read.csv("HiC_loop_counts_5_10_25kb_by_cellline.csv", row.names = 1, check.names = FALSE)

# add CellLine back as a real column
summary_df <- summary_df %>%
  tibble::rownames_to_column("CellLine")

plot_df <- summary_df %>%
  pivot_longer(
    cols = starts_with("loops_"),
    names_to = "Resolution",
    values_to = "LoopCount"
  ) %>%
  mutate(
    Resolution = factor(
      Resolution,
      levels = c("loops_5kb", "loops_10kb", "loops_25kb"),
      labels = c("5 kb", "10 kb", "25 kb")
    ),
    x_num = as.numeric(Resolution)
  ) %>%
  filter(!is.na(LoopCount), LoopCount > 0)  # log10 requires >0

use_xspline <- requireNamespace("ggalt", quietly = TRUE)

p <- ggplot(plot_df, aes(x = x_num, y = LoopCount, color = CellLine, group = CellLine)) +
  (if (use_xspline) ggalt::geom_xspline(spline_shape = 0.6, linewidth = 1.1, alpha = 0.9)
   else geom_line(linewidth = 1.1, alpha = 0.9)) +
  geom_point(size = 2.2, alpha = 0.95) +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("5 kb", "10 kb", "25 kb"),
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  scale_y_log10() +
  labs(
    title = "Hi-C loop counts across resolutions",
    subtitle = "Each line is a cell line; coarser resolution yields fewer detected loops",
    x = "Resolution used for loop calling",
    y = "Number of loops (log10 scale)",
    color = "Cell line"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot",
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ggsave(out_png, p, width = 8.2, height = 5.2, dpi = 300)
ggsave(out_pdf, p, width = 8.2, height = 5.2)
