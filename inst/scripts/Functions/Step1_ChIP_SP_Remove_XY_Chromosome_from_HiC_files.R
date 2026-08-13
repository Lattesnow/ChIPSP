# Remove chrX / chrY rows from ALL *HiC.xls files
library(readxl) 
library(tools)    

hic_files <- list.files(
  path = getwd(),
  pattern = "HiC\\.(xls|xlsx|csv|tsv|txt|bed|tab)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(hic_files) == 0) {
  stop("No *HiC.xls files found in current directory.")
}

for (f in hic_files) {
  
  cat("Processing:", basename(f), "\n")
  
  ext <- tolower(tools::file_ext(f))
  
  hic <- switch(
    ext,
    "csv"  = read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
    "tsv"  = read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
    "txt"  = read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
    "tab"  = read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
    "bed"  = read.delim(f, stringsAsFactors = FALSE, check.names = FALSE, header = FALSE),
    "xlsx" = readxl::read_excel(f) |> as.data.frame(check.names = FALSE),
    "xls"  = tryCatch(
      readxl::read_excel(f) |> as.data.frame(check.names = FALSE),
      error = function(e) {
        read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
      }
    ),
    stop("Unsupported file type: ", ext)
  )
  
  # Detect chromosome columns
  hic_cols <- NULL
  if (all(c("BIN1_CHR", "BIN2_CHR") %in% colnames(hic))) {
    hic_cols <- c("BIN1_CHR", "BIN2_CHR")
  } else if (all(c("BIN1_CHROMOSOME", "BIN2_CHROMOSOME") %in% colnames(hic))) {
    hic_cols <- c("BIN1_CHROMOSOME", "BIN2_CHROMOSOME")
  } else {
    candidates <- c("BIN1_CHR", "BIN2_CHR",
                    "BIN1_CHROMOSOME", "BIN2_CHROMOSOME")
    hic_cols <- intersect(candidates, colnames(hic))
  }
  
  if (length(hic_cols) == 0) {
    warning("  No chromosome columns found — skipped")
    next
  }
  
  xy_regex <- "^(chr)?[XY]$"
  
  keep <- apply(hic[, hic_cols, drop = FALSE], 1, function(x) {
    x <- trimws(as.character(x))
    !any(grepl(xy_regex, x, ignore.case = TRUE) & !is.na(x))
  })
  
  n_removed <- sum(!keep)
  if (n_removed > 0) {
    cat("  Removed", n_removed, "rows containing chrX/chrY\n")
  } else {
    cat("  No chrX/chrY rows found\n")
  }
  
  hic_clean <- hic[keep, , drop = FALSE]
  
  # Output cleaned file as CSV
  out_file <- file.path(
    dirname(f),
    paste0(tools::file_path_sans_ext(basename(f)), "_no_chrX_chrY.csv")
  )
  
  write.csv(
    hic_clean,
    out_file,
    quote = FALSE,
    row.names = FALSE
  )
  
  cat("  Wrote:", basename(out_file), "\n") }
