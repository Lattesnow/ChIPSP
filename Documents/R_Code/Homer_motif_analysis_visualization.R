# ============================================================
#  visualization of HOMER AR motif analysis
# ============================================================

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(Cairo)
library(ggseqlogo)
library(patchwork)

out_dir <- "AR_HOMER_figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1.  theme
# ------------------------------------------------------------

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8, color = "black"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      panel.grid = element_blank()
    )
}

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

count_bed <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- x[x != ""]
  length(x)
}

count_homer_find_hits <- function(file) {
  if (!file.exists(file)) stop("Missing file: ", file)
  
  x <- read.delim(
    file,
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "#",
    check.names = FALSE
  )
  
  if (nrow(x) == 0) return(0)
  
  # HOMER -find output: count unique region IDs from first column
  length(unique(x[[1]]))
}

read_homer_known <- function(file, group_name, top_n = 10) {
  x <- read.delim(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(x) <- make.names(colnames(x))
  
  motif_col <- grep("Motif.Name|Motif", colnames(x), ignore.case = TRUE, value = TRUE)[1]
  p_col <- grep("P.value", colnames(x), ignore.case = TRUE, value = TRUE)[1]
  
  motif_dir <- file.path(dirname(file), "knownResults")
  motif_files <- list.files(motif_dir, pattern = "known[0-9]+\\.motif$", full.names = TRUE)
  motif_files <- motif_files[order(as.numeric(str_extract(basename(motif_files), "[0-9]+")))]
  
  x %>%
    mutate(
      Group = group_name,
      Rank = row_number(),
      Motif_file = motif_files[Rank],
      Motif_raw = .data[[motif_col]],
      Motif = Motif_raw %>%
        str_replace("\\(.*", "") %>%
        str_replace("/.*", "") %>%
        str_trim(),
      P.value = as.numeric(.data[[p_col]]),
      neg_log10_p = -log10(P.value),
      neg_log10_p = ifelse(is.infinite(neg_log10_p), 350, neg_log10_p)
    ) %>%
    arrange(P.value) %>%
    slice_head(n = top_n)
}

read_homer_pwm <- function(motif_file) {
  x <- readLines(motif_file)
  x <- x[!grepl("^>", x)]
  mat <- read.table(text = x, header = FALSE)
  mat <- as.matrix(mat[, 1:4])
  colnames(mat) <- c("A", "C", "G", "T")
  t(mat)
}

plot_homer_motif_dot <- function(df, group_name, out_dir) {
  
  plot_df <- df %>%
    filter(Group == group_name) %>%
    arrange(P.value) %>%
    mutate(
      rank = row_number(),
      Motif_short = str_trunc(Motif, 32),
      Motif_label = paste0(rank, ". ", Motif_short),
      Motif_label = factor(Motif_label, levels = rev(Motif_label))
    )
  
  pwm_list <- lapply(plot_df$Motif_file, read_homer_pwm)
  names(pwm_list) <- plot_df$Motif_label
  
  p_dot <- ggplot(plot_df, aes(x = neg_log10_p, y = Motif_label)) +
    geom_point(aes(size = neg_log10_p, color = neg_log10_p), alpha = 0.95) +
    scale_color_gradientn(
      colors = c("#2F8DD8", "#9B8FBF", "#E46A6A"),
      name = expression(-log[10](P))
    ) +
    scale_size_continuous(range = c(2.5, 7), name = expression(-log[10](P))) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(
      title = paste0(group_name, " motif enrichment"),
      x = expression(-log[10](P)),
      y = NULL
    ) +
    theme_nature() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  p_logo <- ggseqlogo(pwm_list, ncol = 1) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      plot.margin = margin(22, 5, 5, 5)
    )
  
  p_final <- p_dot + p_logo + plot_layout(widths = c(1.25, 1))
  
  safe_group <- gsub("[^A-Za-z0-9]+", "_", group_name)
  
  ggsave(
    file.path(out_dir, paste0("Nature_HOMER_motif_with_logo_", safe_group, ".pdf")),
    p_final,
    device = cairo_pdf,
    width = 10,
    height = 5.5,
    units = "in"
  )
  
  ggsave(
    file.path(out_dir, paste0("Nature_HOMER_motif_with_logo_", safe_group, ".png")),
    p_final,
    width = 10,
    height = 5.5,
    dpi = 600,
    units = "in"
  )
  
  return(p_final)
}
# ------------------------------------------------------------
# 3. AR motif occurrence in three region sets
# ------------------------------------------------------------

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
    Total_regions = count_bed(bed_file),
    AR_motif_regions = count_homer_find_hits(hit_file),
    AR_motif_percent = 100 * AR_motif_regions / Total_regions
  ) %>%
  ungroup() %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "AR original sites",
        "ChIP-SP positive bins",
        "ChIP-SP negative bins"
      )
    )
  )

write_csv(
  motif_occurrence,
  file.path(out_dir, "AR_motif_occurrence_summary.csv")
)

p_occ <- ggplot(motif_occurrence, aes(x = Group, y = AR_motif_percent)) +
  geom_col(width = 0.6, fill = "#3A3A3A") +
  geom_text(
    aes(label = paste0(round(AR_motif_percent, 1), "%")),
    vjust = -0.35,
    size = 2.7
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Canonical AR motif occurrence",
    x = NULL,
    y = "Regions containing AR motif (%)"
  ) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  file.path(out_dir, "Nature_AR_motif_occurrence.pdf"),
  p_occ,
  device = cairo_pdf,
  width = 4.1,
  height = 3.1,
  units = "in"
)

ggsave(
  file.path(out_dir, "Nature_AR_motif_occurrence.png"),
  p_occ,
  width = 4.1,
  height = 3.1,
  dpi = 600,
  units = "in"
)

# ------------------------------------------------------------
# 4. HOMER motif enrichment results
# ------------------------------------------------------------
motif_original <- read_homer_known("homer_AR_original/knownResults.txt", "AR original sites", top_n = 10)
motif_positive <- read_homer_known("homer_HiC_AR_positive/knownResults.txt", "ChIP-SP positive bins", top_n = 10)
motif_negative <- read_homer_known("homer_HiC_AR_negative_random537/knownResults.txt", "ChIP-SP negative bins", top_n = 10)
motif_pos_vs_neg <- read_homer_known("homer_HiC_AR_positive_vs_negative/knownResults.txt", "Positive vs negative", top_n = 10)

motif_all <- bind_rows(motif_original, motif_positive, motif_negative, motif_pos_vs_neg)

plot_homer_motif_dot(motif_all, "AR original sites", out_dir)
plot_homer_motif_dot(motif_all, "ChIP-SP positive bins", out_dir)
plot_homer_motif_dot(motif_all, "ChIP-SP negative bins", out_dir)
plot_homer_motif_dot(motif_all, "Positive vs negative", out_dir)

# ------------------------------------------------------------
# 5. Separate dotplots for each enrichment result
# ------------------------------------------------------------

plot_homer_motif_dot <- function(df, group_name, out_dir) {
  
  plot_df <- df %>%
    filter(Group == group_name) %>%
    arrange(P.value) %>%
    mutate(
      rank = row_number(),
      Motif_short = str_trunc(Motif, 35),
      Motif_label = paste0(rank, ". ", Motif_short),
      Motif_label = factor(Motif_label, levels = rev(Motif_label))
    )
  
  if (nrow(plot_df) == 0) return(NULL)
  
  p <- ggplot(plot_df, aes(x = neg_log10_p, y = Motif_label)) +
    geom_point(
      aes(size = neg_log10_p, color = neg_log10_p),
      alpha = 0.95
    ) +
    scale_color_gradientn(
      colors = c("#2F8DD8", "#9B8FBF", "#E46A6A"),
      name = expression(-log[10](P))
    ) +
    scale_size_continuous(
      range = c(2.5, 7),
      name = expression(-log[10](P))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(group_name, " motif enrichment"),
      x = expression(-log[10](P)),
      y = NULL
    ) +
    theme_nature() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.margin = margin(6, 35, 6, 6),
      legend.position = "right"
    )
  
  safe_group <- gsub("[^A-Za-z0-9]+", "_", group_name)
  
  ggsave(
    file.path(out_dir, paste0("Nature_HOMER_motif_", safe_group, ".pdf")),
    p,
    device = cairo_pdf,
    width = 6.8,
    height = 4.2,
    units = "in"
  )
  
  ggsave(
    file.path(out_dir, paste0("Nature_HOMER_motif_", safe_group, ".png")),
    p,
    width = 6.8,
    height = 4.2,
    dpi = 600,
    units = "in"
  )
  
  return(p)
}

plot_homer_motif_dot(motif_all, "AR original sites", out_dir)
plot_homer_motif_dot(motif_all, "ChIP-SP positive bins", out_dir)
plot_homer_motif_dot(motif_all, "ChIP-SP negative bins", out_dir)
plot_homer_motif_dot(motif_all, "Positive vs negative", out_dir)

# ------------------------------------------------------------
# 6. Print summaries
# ------------------------------------------------------------
motif_all_negative7 <- motif_all %>%
  filter(Group != "ChIP-SP negative bins") %>%
  bind_rows(
    motif_all %>%
      filter(Group == "ChIP-SP negative bins") %>%
      arrange(P.value) %>%
      slice_head(n = 7)
  )

plot_homer_motif_dot(motif_all_negative7, "ChIP-SP negative bins", out_dir)


print(motif_occurrence)

motif_all %>%
  select(Group, Motif_raw, P.value, neg_log10_p) %>%
  print(n = 50)