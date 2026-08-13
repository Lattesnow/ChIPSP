# ============================================================
# Generate absolute negative Hi-C rows for AR ChIP-SP
# ============================================================

library(readxl)
library(dplyr)
library(stringr)
library(readr)

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

hic_file <- "Combined_HiC.xls"
pos_file <- "AR_ChIPSP_mirror_bins.xls"

hic <- read.delim(hic_file, header = TRUE, stringsAsFactors = FALSE)
pos <- read.delim(pos_file, header = TRUE, stringsAsFactors = FALSE)

colnames(hic)
colnames(pos)

# ------------------------------------------------------------
# 2. Define column names
# ------------------------------------------------------------
# Modify these to match your real file column names.

# Hi-C file has two bins per row
hic_chr1_col   <- "BIN1_CHR"
hic_start1_col <- "BIN1_START"
hic_end1_col   <- "BIN1_END"

hic_chr2_col   <- "BIN2_CHR"
hic_start2_col <- "BIN2_START"
hic_end2_col   <- "BIN2_END"

# AR ChIP-SP mirror bin file
pos_chr_col   <- "chr"
pos_start_col <- "start"
pos_end_col   <- "end"

# ------------------------------------------------------------
# 3. Helper function
# ------------------------------------------------------------

clean_chr <- function(x) {
  x <- as.character(x)
  x <- ifelse(str_detect(x, "^chr"), x, paste0("chr", x))
  x
}

make_bin_id <- function(chr, start, end) {
  paste(chr, start, end, sep = "_")
}

# ------------------------------------------------------------
# 4. Standardize AR ChIP-SP positive bins
# ------------------------------------------------------------

pos_bins <- pos %>%
  transmute(
    chr   = clean_chr(.data[[pos_chr_col]]),
    start = as.integer(.data[[pos_start_col]]),
    end   = as.integer(.data[[pos_end_col]])
  ) %>%
  filter(
    !is.na(chr),
    !is.na(start),
    !is.na(end),
    end > start
  ) %>%
  filter(str_detect(chr, "^chr([1-9]|1[0-9]|2[0-2]|X|Y)$")) %>%
  distinct(chr, start, end) %>%
  mutate(pos_bin_id = make_bin_id(chr, start, end))

positive_ids <- pos_bins$pos_bin_id

# ------------------------------------------------------------
# 5. Standardize Hi-C file and create bin IDs for both anchors
# ------------------------------------------------------------

hic2 <- hic %>%
  mutate(
    chr1_clean   = clean_chr(.data[[hic_chr1_col]]),
    start1_clean = as.integer(.data[[hic_start1_col]]),
    end1_clean   = as.integer(.data[[hic_end1_col]]),
    
    chr2_clean   = clean_chr(.data[[hic_chr2_col]]),
    start2_clean = as.integer(.data[[hic_start2_col]]),
    end2_clean   = as.integer(.data[[hic_end2_col]])
  ) %>%
  filter(
    !is.na(chr1_clean), !is.na(start1_clean), !is.na(end1_clean),
    !is.na(chr2_clean), !is.na(start2_clean), !is.na(end2_clean),
    end1_clean > start1_clean,
    end2_clean > start2_clean
  ) %>%
  mutate(
    bin1_id = make_bin_id(chr1_clean, start1_clean, end1_clean),
    bin2_id = make_bin_id(chr2_clean, start2_clean, end2_clean),
    
    bin1_is_AR_ChIPSP_positive = bin1_id %in% positive_ids,
    bin2_is_AR_ChIPSP_positive = bin2_id %in% positive_ids,
    row_has_AR_ChIPSP_positive = bin1_is_AR_ChIPSP_positive | bin2_is_AR_ChIPSP_positive
  )

# ------------------------------------------------------------
# 6. Remove all Hi-C rows containing AR ChIP-SP positive bins
# ------------------------------------------------------------

hic_AR_positive_rows <- hic2 %>%
  filter(row_has_AR_ChIPSP_positive)

hic_AR_negative_rows <- hic2 %>%
  filter(!row_has_AR_ChIPSP_positive)

# ------------------------------------------------------------
# 7. Extract negative bins from remaining negative Hi-C rows
# ------------------------------------------------------------
# These are bins from Hi-C interactions that do not contain any AR ChIP-SP-positive bin.

negative_bin1 <- hic_AR_negative_rows %>%
  transmute(
    chr = chr1_clean,
    start = start1_clean,
    end = end1_clean
  )

negative_bin2 <- hic_AR_negative_rows %>%
  transmute(
    chr = chr2_clean,
    start = start2_clean,
    end = end2_clean
  )

AR_ChIPSP_negative_bins <- bind_rows(negative_bin1, negative_bin2) %>%
  distinct(chr, start, end) %>%
  filter(str_detect(chr, "^chr([1-9]|1[0-9]|2[0-2]|X|Y)$")) %>%
  arrange(chr, start, end) %>%
  mutate(name = paste0("AR_ChIPSP_negative_", row_number()))

AR_ChIPSP_positive_bins <- pos_bins %>%
  select(chr, start, end) %>%
  arrange(chr, start, end) %>%
  mutate(name = paste0("AR_ChIPSP_positive_", row_number()))

# ------------------------------------------------------------
# 8. Save output
# ------------------------------------------------------------

write_tsv(
  AR_ChIPSP_positive_bins %>% select(chr, start, end, name),
  "AR_ChIPSP_positive_bins.bed",
  col_names = FALSE
)

write_tsv(
  AR_ChIPSP_negative_bins %>% select(chr, start, end, name),
  "AR_ChIPSP_negative_bins_absolute.bed",
  col_names = FALSE
)

write_tsv(
  hic_AR_negative_rows,
  "AR_ChIPSP_negative_HiC_rows.tsv"
)

write_tsv(
  hic_AR_positive_rows,
  "AR_ChIPSP_positive_HiC_rows.tsv"
)

# ------------------------------------------------------------
# 9. Filter only 5 kb bins for motif discovery
# ------------------------------------------------------------

AR_ChIPSP_positive_bins_5kb <- AR_ChIPSP_positive_bins %>%
  mutate(width = end - start) %>%
  filter(width == 5000) %>%
  select(chr, start, end, name)

AR_ChIPSP_negative_bins_5kb <- AR_ChIPSP_negative_bins %>%
  mutate(width = end - start) %>%
  filter(width == 5000) %>%
  select(chr, start, end, name)

write_tsv(
  AR_ChIPSP_positive_bins_5kb,
  "AR_ChIPSP_positive_bins_5kb.bed",
  col_names = FALSE
)

write_tsv(
  AR_ChIPSP_negative_bins_5kb,
  "AR_ChIPSP_negative_bins_absolute_5kb.bed",
  col_names = FALSE
)

cat("5 kb AR ChIP-SP positive bins:", nrow(AR_ChIPSP_positive_bins_5kb), "\n")
cat("5 kb AR ChIP-SP negative bins:", nrow(AR_ChIPSP_negative_bins_5kb), "\n")