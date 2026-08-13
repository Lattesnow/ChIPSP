suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# -----------------------------
# Settings
# -----------------------------
base_dir <- "."                 # directory that contains folders like LNCaP/, A549/, ...
hic_subdir <- "HiC"             # inside each cell line folder
res_files <- c("5kb_HiC.xls", "10kb_HiC.xls", "25kb_HiC.xls")

out_csv  <- "HiC_loop_counts_5_10_25kb_by_cellline.csv"
out_png  <- "HiC_loop_counts_5_10_25kb_by_cellline.png"
out_pdf  <- "HiC_loop_counts_5_10_25kb_by_cellline.pdf"

# read hic files
# -----------------------------
read_hic_table <- function(path) {
  # 1) Try readxl if available + file is real Excel
  if (requireNamespace("readxl", quietly = TRUE)) {
    x <- tryCatch(
      readxl::read_excel(path),
      error = function(e) NULL
    )
    if (!is.null(x)) return(as.data.frame(x))
  }
  
  # 2) Fallback: treat as delimited text (common for .xls in pipelines)
  # Try tab first, then comma
  x <- tryCatch(
    read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (!is.null(x)) return(x)
  
  x <- tryCatch(
    read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (!is.null(x)) return(x)
  
  stop("Could not read file as Excel or delimited text: ", path)
}

# Count loops in a HiC file
count_loops <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  df <- read_hic_table(path)
  
  # Most loop files have 1 header row and then data rows.
  # nrow(df) is the correct number of loops (no need -1) as long as header=TRUE.
  nrow(df)
}

cell_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# Keep only those that have HiC subfolder OR contain HiC files
cell_dirs <- cell_dirs[
  sapply(cell_dirs, function(d) dir.exists(file.path(d, hic_subdir)))
]

if (length(cell_dirs) == 0) {
  stop("No cell line folders found that contain a '", hic_subdir, "' subfolder under: ", normalizePath(base_dir))
}

# Build summary table
summary_df <- lapply(cell_dirs, function(cd) {
  cell_line <- basename(cd)
  hic_dir <- file.path(cd, hic_subdir)
  
  paths <- file.path(hic_dir, res_files)
  names(paths) <- res_files
  
  tibble(
    CellLine = cell_line,
    loops_5kb  = count_loops(paths["5kb_HiC.xls"]),
    loops_10kb = count_loops(paths["10kb_HiC.xls"]),
    loops_25kb = count_loops(paths["25kb_HiC.xls"])
  )
}) %>% bind_rows()

# Write CSV with rownames = CellLine (as you requested)
summary_out <- as.data.frame(summary_df)
rownames(summary_out) <- summary_out$CellLine
summary_out$CellLine <- NULL
write.csv(summary_out, out_csv, quote = TRUE)

message("Saved summary CSV: ", out_csv)
