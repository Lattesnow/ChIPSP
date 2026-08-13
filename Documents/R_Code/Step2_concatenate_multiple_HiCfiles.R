# ============================================================
# Merge Hi-C and ChIP files
# Supports CSV, TSV, TXT, TAB, BED, XLS, and XLSX
# Normalizes chromosome names and column names
# Does not delete original files by default
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(tools)
})

wd <- getwd()

# ------------------------------------------------------------
# Universal table reader
# ------------------------------------------------------------

read_any_table <- function(file, file_type = c("HiC", "ChIP")) {
  
  file_type <- match.arg(file_type)
  ext <- tolower(tools::file_ext(file))
  
  message("Reading: ", basename(file))
  
  out <- switch(
    ext,
    
    "csv" = read.csv(
      file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    
    "tsv" = read.delim(
      file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    
    "txt" = read.delim(
      file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    
    "tab" = read.delim(
      file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    
    "bed" = {
      x <- read.delim(
        file,
        header = FALSE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      if (ncol(x) < 3) {
        stop("BED file must contain at least three columns: ", basename(file))
      }
      
      colnames(x)[1:3] <- c("chr", "start", "end")
      
      if (ncol(x) >= 4) {
        colnames(x)[4] <- "pileup"
      }
      
      x
    },
    
    "xls" = tryCatch(
      {
        readxl::read_excel(file) |>
          as.data.frame(check.names = FALSE)
      },
      error = function(e) {
        message(
          "Excel reading failed; attempting tab-delimited reading: ",
          basename(file)
        )
        
        read.delim(
          file,
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    ),
    
    "xlsx" = readxl::read_excel(file) |>
      as.data.frame(check.names = FALSE),
    
    stop("Unsupported file format: ", ext)
  )
  
  out
}

# ------------------------------------------------------------
# Add chr prefix
# ------------------------------------------------------------

add_chr_prefix <- function(x) {
  
  x <- trimws(as.character(x))
  
  need_fix <- (
    !is.na(x) &
      x != "" &
      !grepl("^chr", x, ignore.case = TRUE)
  )
  
  x[need_fix] <- paste0("chr", x[need_fix])
  
  x
}

# ------------------------------------------------------------
# Merge Hi-C files
# ------------------------------------------------------------

hic_files <- list.files(
  path = wd,
  pattern = "HiC.*\\.(xls|xlsx|csv|tsv|txt|tab)$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Exclude already generated combined output
hic_files <- hic_files[
  basename(hic_files) != "Combined_HiC.csv"
]

if (length(hic_files) == 0) {
  
  warning("No Hi-C files found.")
  
} else {
  
  message("Hi-C files detected:")
  
  for (f in hic_files) {
    message("  ", basename(f))
  }
  
  hic_list <- lapply(
    hic_files,
    read_any_table,
    file_type = "HiC"
  )
  
  # Normalize chromosome column names before merging
  hic_list <- lapply(hic_list, function(hic) {
    
    if (
      "BIN1_CHROMOSOME" %in% colnames(hic) &&
      !"BIN1_CHR" %in% colnames(hic)
    ) {
      colnames(hic)[
        colnames(hic) == "BIN1_CHROMOSOME"
      ] <- "BIN1_CHR"
    }
    
    if (
      "BIN2_CHROMOSOME" %in% colnames(hic) &&
      !"BIN2_CHR" %in% colnames(hic)
    ) {
      colnames(hic)[
        colnames(hic) == "BIN2_CHROMOSOME"
      ] <- "BIN2_CHR"
    }
    
    hic
  })
  
  combined_hic <- bind_rows(hic_list)
  
  required_hic_cols <- c(
    "BIN1_CHR",
    "BIN1_START",
    "BIN1_END",
    "BIN2_CHR",
    "BIN2_START",
    "BIN2_END"
  )
  
  missing_hic_cols <- setdiff(
    required_hic_cols,
    colnames(combined_hic)
  )
  
  if (length(missing_hic_cols) > 0) {
    warning(
      "Combined Hi-C data are missing columns: ",
      paste(missing_hic_cols, collapse = ", ")
    )
  }
  
  for (col in intersect(
    c("BIN1_CHR", "BIN2_CHR"),
    colnames(combined_hic)
  )) {
    
    old_values <- combined_hic[[col]]
    
    combined_hic[[col]] <- add_chr_prefix(
      combined_hic[[col]]
    )
    
    n_fixed <- sum(
      old_values != combined_hic[[col]],
      na.rm = TRUE
    )
    
    message(
      "[Hi-C] Added chr prefix in ",
      col,
      " for ",
      n_fixed,
      " rows."
    )
  }
  
  write.csv(
    combined_hic,
    file.path(wd, "Combined_HiC.csv"),
    quote = FALSE,
    row.names = FALSE
  )
  
  message(
    "Created Combined_HiC.csv with ",
    nrow(combined_hic),
    " rows."
  )
}

# ------------------------------------------------------------
# Merge ChIP files
# ------------------------------------------------------------

chip_files <- list.files(
  path = wd,
  pattern = "ChIP.*\\.(bed|xls|xlsx|csv|tsv|txt|tab)$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Exclude already generated combined output
chip_files <- chip_files[
  basename(chip_files) != "Combined_ChIP.csv"
]

if (length(chip_files) == 0) {
  
  warning("No ChIP files found.")
  
} else {
  
  message("ChIP files detected:")
  
  for (f in chip_files) {
    message("  ", basename(f))
  }
  
  chip_list <- lapply(
    chip_files,
    read_any_table,
    file_type = "ChIP"
  )
  
  combined_chip <- bind_rows(chip_list)
  
  # Normalize common alternative column names
  rename_map <- c(
    "chrom" = "chr",
    "chromosome" = "chr",
    "CHR" = "chr",
    "START" = "start",
    "END" = "end",
    "Pileup" = "pileup",
    "PILEUP" = "pileup"
  )
  
  for (old_name in names(rename_map)) {
    
    new_name <- rename_map[[old_name]]
    
    if (
      old_name %in% colnames(combined_chip) &&
      !new_name %in% colnames(combined_chip)
    ) {
      colnames(combined_chip)[
        colnames(combined_chip) == old_name
      ] <- new_name
    }
  }
  
  required_chip_cols <- c(
    "chr",
    "start",
    "end"
  )
  
  missing_chip_cols <- setdiff(
    required_chip_cols,
    colnames(combined_chip)
  )
  
  if (length(missing_chip_cols) > 0) {
    warning(
      "Combined ChIP data are missing columns: ",
      paste(missing_chip_cols, collapse = ", ")
    )
  }
  
  if ("chr" %in% colnames(combined_chip)) {
    
    old_values <- combined_chip$chr
    
    combined_chip$chr <- add_chr_prefix(
      combined_chip$chr
    )
    
    n_fixed <- sum(
      old_values != combined_chip$chr,
      na.rm = TRUE
    )
    
    message(
      "[ChIP] Added chr prefix for ",
      n_fixed,
      " rows."
    )
  }
  
  if (!"pileup" %in% colnames(combined_chip)) {
    combined_chip$pileup <- 1
    
    message(
      "[ChIP] No pileup column detected; assigned pileup = 1."
    )
  }
  
  write.csv(
    combined_chip,
    file.path(wd, "Combined_ChIP.csv"),
    quote = FALSE,
    row.names = FALSE
  )
  
  message(
    "Created Combined_ChIP.csv with ",
    nrow(combined_chip),
    " rows."
  )
}

# ------------------------------------------------------------
# Confirm outputs
# ------------------------------------------------------------

cat("\nOutput check:\n")

cat(
  "Combined_HiC.csv exists:",
  file.exists(file.path(wd, "Combined_HiC.csv")),
  "\n"
)

cat(
  "Combined_ChIP.csv exists:",
  file.exists(file.path(wd, "Combined_ChIP.csv")),
  "\n"
)