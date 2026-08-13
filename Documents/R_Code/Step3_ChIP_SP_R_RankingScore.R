suppressPackageStartupMessages({
  library(data.table)
  library(tools)
})

folder_path <- getwd()

# ============================================================
# 0. Clean intermediate ChIP/Hi-C files
# ============================================================

KEEP_FILES <- c(
  "Combined_ChIP.csv",
  "Combined_HiC.csv"
)

# Only consider genomic input/intermediate files for deletion.
# R scripts, README files, figures, package files, etc. are preserved.
candidate_files <- list.files(
  path = folder_path,
  pattern = "\\.(bed|csv|tsv|txt|tab|xls|xlsx)$",
  full.names = TRUE,
  ignore.case = TRUE
)

candidate_names <- basename(candidate_files)

# Delete only files whose names contain ChIP or HiC,
# excluding the two combined inputs.
files_to_delete <- candidate_files[
  grepl("ChIP|HiC", candidate_names, ignore.case = TRUE) &
    !(tolower(candidate_names) %in% tolower(KEEP_FILES))
]

if (length(files_to_delete) > 0) {
  
  message("The following intermediate ChIP/Hi-C files will be deleted:")
  
  for (f in files_to_delete) {
    message("  ", basename(f))
  }
  
  deletion_status <- file.remove(files_to_delete)
  
  if (any(!deletion_status)) {
    warning(
      "Failed to delete: ",
      paste(
        basename(files_to_delete[!deletion_status]),
        collapse = ", "
      )
    )
  }
  
  message(
    "Deleted ",
    sum(deletion_status),
    " intermediate ChIP/Hi-C files."
  )
  
} else {
  
  message("No intermediate ChIP/Hi-C files required deletion.")
}

# Confirm that required combined files remain
required_files <- file.path(folder_path, KEEP_FILES)

missing_required <- required_files[!file.exists(required_files)]

if (length(missing_required) > 0) {
  stop(
    "Required combined input files are missing:\n",
    paste(basename(missing_required), collapse = "\n")
  )
}

# ============================================================
# 1. Input detection
# ============================================================

chip_path <- file.path(folder_path, "Combined_ChIP.csv")
hic_path  <- file.path(folder_path, "Combined_HiC.csv")

message("ChIP: ", basename(chip_path))
message("HiC : ", basename(hic_path))

FDR_CUTOFF <- 0.05
OVERLAP_MODE <- "any"

# ============================================================
# 2. Universal reader
# ============================================================

read_any_table <- function(file, type = c("ChIP", "HiC")) {
  
  type <- match.arg(type)
  ext <- tolower(tools::file_ext(file))
  
  if (ext == "csv") {
    
    x <- data.table::fread(
      file,
      sep = ",",
      header = TRUE,
      fill = TRUE,
      data.table = TRUE
    )
    
  } else if (ext %in% c("tsv", "txt", "tab")) {
    
    x <- data.table::fread(
      file,
      sep = "\t",
      header = TRUE,
      fill = TRUE,
      data.table = TRUE
    )
    
  } else if (ext == "bed") {
    
    x <- data.table::fread(
      file,
      sep = "\t",
      header = FALSE,
      fill = TRUE,
      data.table = TRUE
    )
    
    if (ncol(x) < 3) {
      stop("BED file must have at least three columns.")
    }
    
    data.table::setnames(
      x,
      old = names(x)[1:3],
      new = c("chr", "start", "end")
    )
    
    if (ncol(x) >= 4 && !"pileup" %in% colnames(x)) {
      data.table::setnames(
        x,
        old = names(x)[4],
        new = "pileup"
      )
    }
    
    if (!"pileup" %in% colnames(x)) {
      x[, pileup := 1]
    }
    
  } else if (ext %in% c("xls", "xlsx")) {
    
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(
        "Package `readxl` is required to read Excel files."
      )
    }
    
    x <- data.table::as.data.table(
      readxl::read_excel(file)
    )
    
  } else {
    
    stop("Unsupported file type: ", ext)
  }
  
  x
}

# ============================================================
# 3. Read tables
# ============================================================

chip <- read_any_table(
  chip_path,
  type = "ChIP"
)

hic <- read_any_table(
  hic_path,
  type = "HiC"
)

message(
  "Loaded ChIP: ",
  nrow(chip),
  " rows and ",
  ncol(chip),
  " columns."
)

message(
  "Loaded Hi-C: ",
  nrow(hic),
  " rows and ",
  ncol(hic),
  " columns."
)

# ============================================================
# 4. Normalize ChIP columns
# ============================================================

# Normalize common alternative column names
chip_name_map <- c(
  "chrom" = "chr",
  "chromosome" = "chr",
  "CHR" = "chr",
  "START" = "start",
  "END" = "end",
  "Pileup" = "pileup",
  "PILEUP" = "pileup"
)

for (old_name in names(chip_name_map)) {
  
  new_name <- chip_name_map[[old_name]]
  
  if (
    old_name %in% colnames(chip) &&
    !new_name %in% colnames(chip)
  ) {
    data.table::setnames(
      chip,
      old = old_name,
      new = new_name
    )
  }
}

required_chip_cols <- c(
  "chr",
  "start",
  "end"
)

missing_chip_cols <- setdiff(
  required_chip_cols,
  colnames(chip)
)

if (length(missing_chip_cols) > 0) {
  stop(
    "ChIP file is missing required columns: ",
    paste(missing_chip_cols, collapse = ", ")
  )
}

if (!"pileup" %in% colnames(chip)) {
  chip[, pileup := 1]
  
  message(
    "No `pileup` column found; assigned pileup = 1 to all peaks."
  )
}

chip <- chip[, .(
  chr,
  start,
  end,
  pileup
)]

# ============================================================
# 5. Normalize Hi-C columns
# ============================================================

if (
  "BIN1_CHROMOSOME" %in% colnames(hic) &&
  !"BIN1_CHR" %in% colnames(hic)
) {
  data.table::setnames(
    hic,
    "BIN1_CHROMOSOME",
    "BIN1_CHR"
  )
}

if (
  "BIN2_CHROMOSOME" %in% colnames(hic) &&
  !"BIN2_CHR" %in% colnames(hic)
) {
  data.table::setnames(
    hic,
    "BIN2_CHROMOSOME",
    "BIN2_CHR"
  )
}

required_hic_cols <- c(
  "BIN1_CHR",
  "BIN1_START",
  "BIN1_END",
  "BIN2_CHR",
  "BIN2_START",
  "BIN2_END",
  "FDR"
)

missing_hic_cols <- setdiff(
  required_hic_cols,
  colnames(hic)
)

if (length(missing_hic_cols) > 0) {
  stop(
    "Hi-C file is missing required columns: ",
    paste(missing_hic_cols, collapse = ", ")
  )
}

hic <- hic[, .(
  BIN1_CHR,
  BIN1_START,
  BIN1_END,
  BIN2_CHR,
  BIN2_START,
  BIN2_END,
  FDR
)]

# ============================================================
# 6. Normalize chromosome names
# ============================================================

normalize_chr <- function(x) {
  
  x <- trimws(as.character(x))
  
  needs_prefix <- (
    !is.na(x) &
      x != "" &
      !grepl("^chr", x, ignore.case = TRUE)
  )
  
  x[needs_prefix] <- paste0(
    "chr",
    x[needs_prefix]
  )
  
  x <- sub(
    "^CHR",
    "chr",
    x,
    ignore.case = TRUE
  )
  
  x
}

chip[, chr := normalize_chr(chr)]

hic[, `:=`(
  BIN1_CHR = normalize_chr(BIN1_CHR),
  BIN2_CHR = normalize_chr(BIN2_CHR)
)]

# ============================================================
# 7. Type conversion
# ============================================================

chip[, `:=`(
  start = suppressWarnings(as.integer(start)),
  end = suppressWarnings(as.integer(end)),
  pileup = suppressWarnings(as.numeric(pileup))
)]

hic[, `:=`(
  BIN1_START = suppressWarnings(as.integer(BIN1_START)),
  BIN1_END = suppressWarnings(as.integer(BIN1_END)),
  BIN2_START = suppressWarnings(as.integer(BIN2_START)),
  BIN2_END = suppressWarnings(as.integer(BIN2_END)),
  FDR = suppressWarnings(as.numeric(FDR))
)]

# Remove invalid ChIP rows
valid_chip <- (
  !is.na(chip$chr) &
    chip$chr != "" &
    !is.na(chip$start) &
    !is.na(chip$end) &
    chip$end > chip$start &
    !is.na(chip$pileup)
)

n_invalid_chip <- sum(!valid_chip)

if (n_invalid_chip > 0) {
  warning(
    "Removed ",
    n_invalid_chip,
    " invalid ChIP rows."
  )
}

chip <- chip[valid_chip]

# Remove invalid Hi-C rows
valid_hic <- (
  !is.na(hic$BIN1_CHR) &
    hic$BIN1_CHR != "" &
    !is.na(hic$BIN1_START) &
    !is.na(hic$BIN1_END) &
    hic$BIN1_END > hic$BIN1_START &
    !is.na(hic$BIN2_CHR) &
    hic$BIN2_CHR != "" &
    !is.na(hic$BIN2_START) &
    !is.na(hic$BIN2_END) &
    hic$BIN2_END > hic$BIN2_START &
    !is.na(hic$FDR)
)

n_invalid_hic <- sum(!valid_hic)

if (n_invalid_hic > 0) {
  warning(
    "Removed ",
    n_invalid_hic,
    " invalid Hi-C rows."
  )
}

hic <- hic[valid_hic]

message(
  "Rows: ChIP = ",
  nrow(chip),
  "; Hi-C before filtering = ",
  nrow(hic)
)

# ============================================================
# 8. Optional Hi-C FDR filter
# ============================================================

if (
  !is.null(FDR_CUTOFF) &&
  length(FDR_CUTOFF) == 1 &&
  !is.na(FDR_CUTOFF) &&
  is.finite(FDR_CUTOFF) &&
  FDR_CUTOFF < 1
) {
  hic <- hic[FDR <= FDR_CUTOFF]
}

message(
  "Rows: Hi-C after FDR filtering = ",
  nrow(hic)
)

if (nrow(hic) == 0) {
  stop(
    "No Hi-C loops remain after FDR filtering."
  )
}

# ============================================================
# 9. Overlap ChIP with Hi-C BIN1 and project to BIN2
# ============================================================

data.table::setkey(
  chip,
  chr,
  start,
  end
)

bin1 <- hic[, .(
  chr = BIN1_CHR,
  start = BIN1_START,
  end = BIN1_END,
  partner_chr = BIN2_CHR,
  partner_start = BIN2_START,
  partner_end = BIN2_END,
  FDR
)]

data.table::setkey(
  bin1,
  chr,
  start,
  end
)

ov1 <- data.table::foverlaps(
  bin1,
  chip,
  type = OVERLAP_MODE,
  nomatch = 0L
)

new_df1 <- ov1[, .(
  chr = partner_chr,
  start = partner_start,
  end = partner_end,
  pileup,
  FDR,
  source_anchor = "BIN1"
)]

# ============================================================
# 10. Overlap ChIP with Hi-C BIN2 and project to BIN1
# ============================================================

bin2 <- hic[, .(
  chr = BIN2_CHR,
  start = BIN2_START,
  end = BIN2_END,
  partner_chr = BIN1_CHR,
  partner_start = BIN1_START,
  partner_end = BIN1_END,
  FDR
)]

data.table::setkey(
  bin2,
  chr,
  start,
  end
)

ov2 <- data.table::foverlaps(
  bin2,
  chip,
  type = OVERLAP_MODE,
  nomatch = 0L
)

new_df2 <- ov2[, .(
  chr = partner_chr,
  start = partner_start,
  end = partner_end,
  pileup,
  FDR,
  source_anchor = "BIN2"
)]

message(
  "Projected rows: BIN1 = ",
  nrow(new_df1),
  "; BIN2 = ",
  nrow(new_df2)
)


# ============================================================
# 11. Remove within-process duplicates, then merge
# ============================================================

n_df1_before <- nrow(new_df1)
n_df2_before <- nrow(new_df2)

# Remove duplicate rows generated within the BIN1 process only
new_df1 <- unique(
  new_df1,
  by = c(
    "chr",
    "start",
    "end",
    "pileup",
    "FDR",
    "source_anchor"
  )
)

# Remove duplicate rows generated within the BIN2 process only
new_df2 <- unique(
  new_df2,
  by = c(
    "chr",
    "start",
    "end",
    "pileup",
    "FDR",
    "source_anchor"
  )
)

message(
  "BIN1 process: removed ",
  n_df1_before - nrow(new_df1),
  " duplicated rows."
)

message(
  "BIN2 process: removed ",
  n_df2_before - nrow(new_df2),
  " duplicated rows."
)

# Combine the two processes

final_matrix <- data.table::rbindlist(
  list(
    new_df1,
    new_df2
  ),
  use.names = TRUE,
  fill = TRUE
)

message(
  "Final rows after within-process duplicate removal: ",
  nrow(final_matrix)
)

if (nrow(final_matrix) == 0) {
  stop(
    "No ChIP peaks overlapped any retained Hi-C loop anchors."
  )
}

# ------------------------------------------------------------
# Min-max normalization
# ------------------------------------------------------------

rng01 <- function(x) {
  
  x <- as.numeric(x)
  valid <- x[is.finite(x)]
  
  if (length(valid) == 0) {
    return(rep(NA_real_, length(x)))
  }
  
  r <- range(valid, na.rm = TRUE)
  
  if (isTRUE(all.equal(r[1], r[2]))) {
    return(rep(0, length(x)))
  }
  
  (x - r[1]) / (r[2] - r[1])
}

final_matrix[, pileup_norm := rng01(pileup)]
final_matrix[, fdr_norm := rng01(FDR)]
final_matrix[, score := pileup_norm - fdr_norm]

data.table::setorder(
  final_matrix,
  -score,
  FDR,
  -pileup
)

final_matrix[, rank := seq_len(.N)]

# ============================================================
# 12. Write output
# ============================================================

out_file <- file.path(
  folder_path,
  "final_ranked_output.csv"
)

data.table::fwrite(
  final_matrix,
  file = out_file,
  sep = ",",
  quote = FALSE
)

message(
  "Wrote ",
  basename(out_file),
  " with ",
  nrow(final_matrix),
  " rows."
)