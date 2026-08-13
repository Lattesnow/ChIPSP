# ============================================================
#AR motif occurrence density figure only
# ============================================================

library(dplyr)
library(ggplot2)
library(readr)
library(Cairo)

out_dir <- "AR_motif_occurrence_figure"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8, color = "black"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      panel.grid = element_blank()
    )
}

read_bed_safe <- function(file) {
  bed <- read.delim(
    file,
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "#",
    check.names = FALSE
  )
  
  bed <- bed[, 1:3]
  colnames(bed) <- c("chr", "start", "end")
  
  bed %>%
    mutate(
      start = suppressWarnings(as.numeric(start)),
      end = suppressWarnings(as.numeric(end))
    ) %>%
    filter(!is.na(start), !is.na(end), end > start)
}

count_bed_rows <- function(file) {
  nrow(read_bed_safe(file))
}

count_bed_bp <- function(file) {
  bed <- read_bed_safe(file)
  sum(bed$end - bed$start)
}

count_homer_find_hits <- function(file) {
  x <- read.delim(
    file,
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "#",
    check.names = FALSE
  )
  
  if (nrow(x) == 0) return(0)
  
  # one HOMER -find row = one motif hit
  nrow(x)
}

occurrence_input <- tibble::tribble(
  ~Group, ~bed_file, ~hit_file,
  "AR original sites", 
  "AR_original_sites.bed", 
  "AR_original_ARmotif_hits.txt",
  
  "ChIP-SP positive bins", 
  "AR_ChIPSP_positive_bins.bed", 
  "HiC_AR_ChIPSP_positive_ARmotif_hits.txt",
  
  "ChIP-SP negative bins", 
  "AR_ChIPSP_negative_bins_absolute.bed", 
  "HiC_AR_ChIPSP_negative_ARmotif_hits.txt"
)

motif_occurrence <- occurrence_input %>%
  rowwise() %>%
  mutate(
    Total_regions = count_bed_rows(bed_file),
    Total_bp = count_bed_bp(bed_file),
    AR_motif_hits = count_homer_find_hits(hit_file),
    AR_motif_hits_per_kb = AR_motif_hits / Total_bp * 1000
  ) %>%
  ungroup() %>%
  mutate(
    Group = factor(
      Group,
      levels = c("AR original sites", "ChIP-SP positive bins", "ChIP-SP negative bins")
    )
  )

write_csv(
  motif_occurrence,
  file.path(out_dir, "AR_motif_occurrence_density_summary.csv")
)

p_occ <- ggplot(motif_occurrence, aes(x = Group, y = AR_motif_hits_per_kb)) +
  geom_col(width = 0.6, fill = "#3A3A3A") +
  geom_text(
    aes(label = round(AR_motif_hits_per_kb, 3)),
    vjust = -0.35,
    size = 2.7
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Canonical AR motif density",
    x = NULL,
    y = "AR motif hits per kb"
  ) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  file.path(out_dir, "Nature_AR_motif_density.pdf"),
  p_occ,
  device = cairo_pdf,
  width = 4.2,
  height = 3.1,
  units = "in"
)

ggsave(
  file.path(out_dir, "Nature_AR_motif_density.png"),
  p_occ,
  width = 4.2,
  height = 3.1,
  dpi = 600,
  units = "in"
)

print(motif_occurrence)