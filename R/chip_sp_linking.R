#' ChIP-SP Core Spatial Integration and Ranking
#'
#' @description
#' Integrates ChIP-seq peaks with Hi-C chromatin loops by identifying
#' overlaps between ChIP-seq peaks and either Hi-C loop anchor. When a
#' ChIP-seq peak overlaps one loop anchor, the partner anchor is reported
#' as the spatially linked ChIP-SP region.
#'
#' Spatially linked regions are ranked using normalized ChIP-seq pileup
#' and Hi-C loop false-discovery rate:
#'
#' \deqn{
#' score = normalized\ pileup - normalized\ FDR
#' }
#'
#' Duplicate projected rows are removed independently within the BIN1
#' projection process and within the BIN2 projection process. Identical
#' rows produced independently by BIN1 and BIN2 are retained as separate
#' output rows.
#'
#' The function uses \code{data.table::foverlaps()} for efficient genomic
#' interval matching and avoids construction of a full Cartesian product.
#'
#' @param chip_file Character scalar. Path to a ChIP-seq peak file.
#'   Supported formats are BED, CSV, TSV, TXT, TAB, XLS, and XLSX.
#'   Delimited files must contain \code{chr}, \code{start}, and \code{end}
#'   columns. A \code{pileup} column is optional; when absent, all peaks
#'   are assigned a pileup value of 1.
#'
#' @param hic_df A data.frame or data.table containing Hi-C loops,
#'   typically returned by \code{mergeHiCLoops()}. Required columns are
#'   \code{BIN1_CHR}, \code{BIN1_START}, \code{BIN1_END},
#'   \code{BIN2_CHR}, \code{BIN2_START}, \code{BIN2_END}, and \code{FDR}.
#'
#' @param fdr_cutoff Numeric scalar specifying the maximum Hi-C loop FDR
#'   retained before overlap analysis. Default is \code{0.05}. Set to
#'   \code{1} or \code{NULL} to disable FDR filtering.
#'
#' @param overlap_mode Character scalar passed to
#'   \code{data.table::foverlaps(type = ...)}. Default is \code{"any"}.
#'
#' @param add_chr_prefix Logical. If \code{TRUE}, chromosome values lacking
#'   the \code{"chr"} prefix are normalized before overlap analysis.
#'   Default is \code{TRUE}.
#'
#' @return A data.frame containing spatially linked genomic regions ranked
#'   by ChIP-SP score. Output columns include:
#'
#' \itemize{
#'   \item \code{chr}, \code{start}, and \code{end}: partner-anchor coordinates;
#'   \item \code{pileup}: ChIP-seq peak signal;
#'   \item \code{FDR}: Hi-C loop confidence;
#'   \item \code{source_anchor}: loop anchor overlapped by the ChIP-seq peak;
#'   \item \code{pileup_norm}: min-max normalized pileup;
#'   \item \code{fdr_norm}: min-max normalized FDR;
#'   \item \code{score}: final ChIP-SP ranking score;
#'   \item \code{rank}: descending rank based on ChIP-SP score.
#' }
#'
#' @examples
#' hic_files <- list.files(
#'   pattern = "HiC.*\\.(csv|tsv|txt|tab|xls|xlsx)$",
#'   full.names = TRUE,
#'   ignore.case = TRUE
#' )
#'
#' hic_df <- mergeHiCLoops(hic_files)
#'
#' result <- chipSPLink(
#'   chip_file = "AR_ChIP.bed",
#'   hic_df = hic_df,
#'   fdr_cutoff = 0.05
#' )
#'
#' head(result)
#'
#' @export
chipSPLink <- function(
    chip_file = NULL,
    hic_df,
    chip_path = getwd(),
    fdr_cutoff = 0.05,
    overlap_mode = "any",
    add_chr_prefix = TRUE) {

  # ----------------------------------------------------------
  # 1. Validate inputs
  # ----------------------------------------------------------

  if (!is.null(chip_file)) {

    if (
      !is.character(chip_file) ||
        length(chip_file) != 1L ||
        is.na(chip_file) ||
        chip_file == ""
    ) {
      stop(
        "`chip_file` must be NULL or one valid file path."
      )
    }

    if (!file.exists(chip_file)) {
      stop(
        "ChIP file not found: ",
        chip_file
      )
    }
  }

  if (!is.data.frame(hic_df)) {
    stop(
      "`hic_df` must be a data.frame or data.table, ",
      "typically returned by `mergeHiCLoops()`."
    )
  }

  valid_overlap_modes <- c(
    "any",
    "within",
    "start",
    "end",
    "equal"
  )

  if (
    !is.character(overlap_mode) ||
      length(overlap_mode) != 1L ||
      !overlap_mode %in% valid_overlap_modes
  ) {
    stop(
      "`overlap_mode` must be one of: ",
      paste(valid_overlap_modes, collapse = ", ")
    )
  }

  if (
    !is.logical(add_chr_prefix) ||
      length(add_chr_prefix) != 1L ||
      is.na(add_chr_prefix)
  ) {
    stop("`add_chr_prefix` must be TRUE or FALSE.")
  }

  # ----------------------------------------------------------
  # 2. Read ChIP-seq input
  # ----------------------------------------------------------

  chip <- readChIPFile(
    file = chip_file,
    path = chip_path
  )

  # data.table conversion may otherwise modify the input object by reference
  hic <- data.table::as.data.table(
    data.table::copy(hic_df)
  )

  # ----------------------------------------------------------
  # 3. Normalize ChIP column names
  # ----------------------------------------------------------

  chip_names <- trimws(colnames(chip))
  chip_names_upper <- toupper(chip_names)

  chip_name_map <- c(
    CHR = "chr",
    CHROM = "chr",
    CHROMOSOME = "chr",
    START = "start",
    END = "end",
    PILEUP = "pileup"
  )

  matched_chip_names <- (
    chip_names_upper %in% names(chip_name_map)
  )

  chip_names[matched_chip_names] <- unname(
    chip_name_map[chip_names_upper[matched_chip_names]]
  )

  data.table::setnames(
    chip,
    old = colnames(chip),
    new = chip_names
  )

  duplicated_chip_names <- unique(
    colnames(chip)[duplicated(colnames(chip))]
  )

  if (length(duplicated_chip_names) > 0L) {
    stop(
      "ChIP column-name normalization created duplicate columns: ",
      paste(duplicated_chip_names, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 4. Normalize Hi-C column names
  # ----------------------------------------------------------

  hic_names <- trimws(colnames(hic))
  hic_names_upper <- toupper(hic_names)

  hic_name_map <- c(
    BIN1_CHROMOSOME = "BIN1_CHR",
    BIN1_CHR = "BIN1_CHR",
    BIN1_START = "BIN1_START",
    BIN1_END = "BIN1_END",
    BIN2_CHROMOSOME = "BIN2_CHR",
    BIN2_CHR = "BIN2_CHR",
    BIN2_START = "BIN2_START",
    BIN2_END = "BIN2_END",
    FDR = "FDR"
  )

  matched_hic_names <- (
    hic_names_upper %in% names(hic_name_map)
  )

  hic_names[matched_hic_names] <- unname(
    hic_name_map[hic_names_upper[matched_hic_names]]
  )

  data.table::setnames(
    hic,
    old = colnames(hic),
    new = hic_names
  )

  duplicated_hic_names <- unique(
    colnames(hic)[duplicated(colnames(hic))]
  )

  if (length(duplicated_hic_names) > 0L) {
    stop(
      "Hi-C column-name normalization created duplicate columns: ",
      paste(duplicated_hic_names, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 5. Validate required columns
  # ----------------------------------------------------------

  required_chip_columns <- c(
    "chr",
    "start",
    "end"
  )

  missing_chip_columns <- setdiff(
    required_chip_columns,
    colnames(chip)
  )

  if (length(missing_chip_columns) > 0L) {
    stop(
      "Missing ChIP columns: ",
      paste(missing_chip_columns, collapse = ", ")
    )
  }

  if (!"pileup" %in% colnames(chip)) {
    chip[, pileup := 1]

    message(
      "No `pileup` column detected. ",
      "All ChIP peaks were assigned pileup = 1."
    )
  }

  required_hic_columns <- c(
    "BIN1_CHR",
    "BIN1_START",
    "BIN1_END",
    "BIN2_CHR",
    "BIN2_START",
    "BIN2_END",
    "FDR"
  )

  missing_hic_columns <- setdiff(
    required_hic_columns,
    colnames(hic)
  )

  if (length(missing_hic_columns) > 0L) {
    stop(
      "Missing Hi-C columns: ",
      paste(missing_hic_columns, collapse = ", ")
    )
  }

  chip <- chip[, .(
    chr,
    start,
    end,
    pileup
  )]

  hic <- hic[, .(
    BIN1_CHR,
    BIN1_START,
    BIN1_END,
    BIN2_CHR,
    BIN2_START,
    BIN2_END,
    FDR
  )]

  # ----------------------------------------------------------
  # 6. Convert data types
  # ----------------------------------------------------------

  chip[, `:=`(
    chr = trimws(as.character(chr)),
    start = suppressWarnings(as.integer(start)),
    end = suppressWarnings(as.integer(end)),
    pileup = suppressWarnings(as.numeric(pileup))
  )]

  hic[, `:=`(
    BIN1_CHR = trimws(as.character(BIN1_CHR)),
    BIN1_START = suppressWarnings(as.integer(BIN1_START)),
    BIN1_END = suppressWarnings(as.integer(BIN1_END)),
    BIN2_CHR = trimws(as.character(BIN2_CHR)),
    BIN2_START = suppressWarnings(as.integer(BIN2_START)),
    BIN2_END = suppressWarnings(as.integer(BIN2_END)),
    FDR = suppressWarnings(as.numeric(FDR))
  )]

  # ----------------------------------------------------------
  # 7. Remove invalid rows
  # ----------------------------------------------------------

  n_chip_before_validation <- nrow(chip)

  chip <- chip[
    !is.na(chr) &
      chr != "" &
      !is.na(start) &
      !is.na(end) &
      end > start &
      !is.na(pileup)
  ]

  n_chip_removed <- (
    n_chip_before_validation - nrow(chip)
  )

  if (n_chip_removed > 0L) {
    warning(
      "Removed ",
      n_chip_removed,
      " invalid ChIP rows."
    )
  }

  n_hic_before_validation <- nrow(hic)

  hic <- hic[
    !is.na(BIN1_CHR) &
      BIN1_CHR != "" &
      !is.na(BIN1_START) &
      !is.na(BIN1_END) &
      BIN1_END > BIN1_START &
      !is.na(BIN2_CHR) &
      BIN2_CHR != "" &
      !is.na(BIN2_START) &
      !is.na(BIN2_END) &
      BIN2_END > BIN2_START &
      !is.na(FDR)
  ]

  n_hic_removed <- (
    n_hic_before_validation - nrow(hic)
  )

  if (n_hic_removed > 0L) {
    warning(
      "Removed ",
      n_hic_removed,
      " invalid Hi-C rows."
    )
  }

  # ----------------------------------------------------------
  # 8. Normalize chromosome prefixes
  # ----------------------------------------------------------

  normalize_chr <- function(x) {

    x <- trimws(as.character(x))

    missing_value <- is.na(x) | x == ""

    needs_prefix <- (
      !missing_value &
        !grepl("^chr", x, ignore.case = TRUE)
    )

    x[needs_prefix] <- paste0(
      "chr",
      x[needs_prefix]
    )

    x[!missing_value] <- sub(
      "^chr",
      "chr",
      x[!missing_value],
      ignore.case = TRUE
    )

    x
  }

  if (add_chr_prefix) {

    chip[, chr := normalize_chr(chr)]

    hic[, `:=`(
      BIN1_CHR = normalize_chr(BIN1_CHR),
      BIN2_CHR = normalize_chr(BIN2_CHR)
    )]
  }

  message(
    "Rows: ChIP = ",
    nrow(chip),
    "; Hi-C before FDR filtering = ",
    nrow(hic)
  )

  # ----------------------------------------------------------
  # 9. Filter Hi-C loops by FDR
  # ----------------------------------------------------------

  if (!is.null(fdr_cutoff)) {

    if (
      !is.numeric(fdr_cutoff) ||
        length(fdr_cutoff) != 1L ||
        is.na(fdr_cutoff) ||
        !is.finite(fdr_cutoff) ||
        fdr_cutoff < 0
    ) {
      stop(
        "`fdr_cutoff` must be NULL or one finite ",
        "non-negative numeric value."
      )
    }

    if (fdr_cutoff < 1) {
      hic <- hic[FDR <= fdr_cutoff]
    }
  }

  message(
    "Rows: Hi-C after FDR filtering = ",
    nrow(hic)
  )

  if (nrow(chip) == 0L) {
    stop(
      "No valid ChIP-seq peaks remain after input processing."
    )
  }

  if (nrow(hic) == 0L) {
    warning(
      "No Hi-C loops remain after input processing."
    )

    return(data.frame())
  }

  # ----------------------------------------------------------
  # 10. Overlap ChIP peaks with BIN1 and project to BIN2
  # ----------------------------------------------------------

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

  overlap_bin1 <- data.table::foverlaps(
    bin1,
    chip,
    type = overlap_mode,
    nomatch = 0L
  )

  projected_from_bin1 <- overlap_bin1[, .(
    chr = partner_chr,
    start = partner_start,
    end = partner_end,
    pileup,
    FDR,
    source_anchor = "BIN1"
  )]

  # ----------------------------------------------------------
  # 11. Overlap ChIP peaks with BIN2 and project to BIN1
  # ----------------------------------------------------------

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

  overlap_bin2 <- data.table::foverlaps(
    bin2,
    chip,
    type = overlap_mode,
    nomatch = 0L
  )

  projected_from_bin2 <- overlap_bin2[, .(
    chr = partner_chr,
    start = partner_start,
    end = partner_end,
    pileup,
    FDR,
    source_anchor = "BIN2"
  )]

  message(
    "Projected rows before within-process duplicate removal: ",
    "BIN1 = ",
    nrow(projected_from_bin1),
    "; BIN2 = ",
    nrow(projected_from_bin2)
  )

  # ----------------------------------------------------------
  # 12. Remove duplicates within each projection process
  # ----------------------------------------------------------

  n_bin1_before_deduplication <- nrow(
    projected_from_bin1
  )

  n_bin2_before_deduplication <- nrow(
    projected_from_bin2
  )

  projected_from_bin1 <- unique(
    projected_from_bin1,
    by = c(
      "chr",
      "start",
      "end",
      "pileup",
      "FDR",
      "source_anchor"
    )
  )

  projected_from_bin2 <- unique(
    projected_from_bin2,
    by = c(
      "chr",
      "start",
      "end",
      "pileup",
      "FDR",
      "source_anchor"
    )
  )

  n_bin1_removed <- (
    n_bin1_before_deduplication -
      nrow(projected_from_bin1)
  )

  n_bin2_removed <- (
    n_bin2_before_deduplication -
      nrow(projected_from_bin2)
  )

  message(
    "BIN1 projection: removed ",
    n_bin1_removed,
    " within-process duplicated rows."
  )

  message(
    "BIN2 projection: removed ",
    n_bin2_removed,
    " within-process duplicated rows."
  )

  # ----------------------------------------------------------
  # 13. Merge independently generated projection sets
  # ----------------------------------------------------------

  final_matrix <- data.table::rbindlist(
    list(
      projected_from_bin1,
      projected_from_bin2
    ),
    use.names = TRUE,
    fill = TRUE
  )

  # Do not deduplicate final_matrix.
  #
  # An identical row generated independently through BIN1 and BIN2
  # must remain as two separate rows because it originated from two
  # different projection processes.

  message(
    "Projected rows after within-process duplicate removal: ",
    "BIN1 = ",
    nrow(projected_from_bin1),
    "; BIN2 = ",
    nrow(projected_from_bin2),
    "; combined = ",
    nrow(final_matrix)
  )

  if (nrow(final_matrix) == 0L) {
    warning(
      "No ChIP-seq peaks overlapped retained Hi-C loop anchors."
    )

    return(as.data.frame(final_matrix))
  }

  # ----------------------------------------------------------
  # 14. Calculate ChIP-SP ranking score
  # ----------------------------------------------------------

  normalize_01 <- function(x) {

    x <- as.numeric(x)

    valid_values <- x[is.finite(x)]

    if (length(valid_values) == 0L) {
      return(rep(NA_real_, length(x)))
    }

    value_range <- range(
      valid_values,
      na.rm = TRUE
    )

    if (
      isTRUE(
        all.equal(
          value_range[1],
          value_range[2]
        )
      )
    ) {
      return(rep(0, length(x)))
    }

    (x - value_range[1]) /
      (value_range[2] - value_range[1])
  }

  final_matrix[
    ,
    pileup_norm := normalize_01(pileup)
  ]

  final_matrix[
    ,
    fdr_norm := normalize_01(FDR)
  ]

  final_matrix[
    ,
    score := pileup_norm - fdr_norm
  ]

  data.table::setorder(
    final_matrix,
    -score,
    FDR,
    -pileup
  )

  final_matrix[, rank := seq_len(.N)]

  message(
    "Final ChIP-SP output: ",
    nrow(final_matrix),
    " spatially linked rows."
  )

  as.data.frame(final_matrix)
}
