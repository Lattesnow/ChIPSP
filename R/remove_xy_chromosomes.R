#' Remove chrX/chrY from Hi-C loop files
#'
#' Reads one or more Hi-C loop files, removes interactions involving
#' chromosome X or chromosome Y, and writes cleaned files in CSV format.
#' Supported input formats include CSV, TSV, TXT, TAB, BED, XLS, and XLSX.
#' Chromosome columns are detected automatically.
#'
#' If `files` is not supplied, the function automatically searches for
#' Hi-C files in `path`.
#'
#' @param files Optional character vector of Hi-C file paths.
#'   If NULL, matching Hi-C files are automatically detected in `path`.
#' @param path Directory to search when `files = NULL`.
#'   Default is the current working directory.
#'
#' @return Invisibly returns a character vector containing the output file names.
#' @export
#'
#' @examples
#' # Automatically detect Hi-C files in current directory
#' removeXYChromosomes()
#'
#' # Automatically detect Hi-C files in another directory
#' removeXYChromosomes(path = "HiC_results")
#'
#' # Or provide files manually
#' hic_files <- list.files(
#'   pattern = "HiC\\.(csv|tsv|txt|xls|xlsx)$",
#'   full.names = TRUE
#' )
#'
#' removeXYChromosomes(hic_files)
removeXYChromosomes <- function(files = NULL, path = getwd()) {

  # ----------------------------------------------------------
  # Automatically detect Hi-C files if files are not supplied
  # ----------------------------------------------------------
  if (is.null(files)) {

    files <- list.files(
      path = path,
      pattern = "HiC\\.(xls|xlsx|csv|tsv|txt|bed|tab)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  }

  if (length(files) == 0) {
    stop(
      "No Hi-C files found. ",
      "Expected filenames ending with HiC.xls, HiC.xlsx, HiC.csv, ",
      "HiC.tsv, HiC.txt, HiC.bed, or HiC.tab."
    )
  }

  output_files <- character(length(files))

  for (i in seq_along(files)) {

    f <- files[i]

    message("Processing: ", basename(f))

    ext <- tolower(tools::file_ext(f))

    hic <- switch(
      ext,

      csv = read.csv(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),

      tsv = read.delim(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),

      txt = read.delim(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),

      tab = read.delim(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),

      bed = read.delim(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        header = FALSE
      ),

      xlsx = readxl::read_excel(f) |>
        as.data.frame(
          check.names = FALSE
        ),

      xls = tryCatch(
        readxl::read_excel(f) |>
          as.data.frame(
            check.names = FALSE
          ),

        error = function(e) {
          read.delim(
            f,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        }
      ),

      stop(
        "Unsupported file type: ",
        ext
      )
    )

    # --------------------------------------------------------
    # Detect chromosome columns
    # --------------------------------------------------------
    if (
      all(
        c(
          "BIN1_CHR",
          "BIN2_CHR"
        ) %in% colnames(hic)
      )
    ) {

      hic_cols <- c(
        "BIN1_CHR",
        "BIN2_CHR"
      )

    } else if (
      all(
        c(
          "BIN1_CHROMOSOME",
          "BIN2_CHROMOSOME"
        ) %in% colnames(hic)
      )
    ) {

      hic_cols <- c(
        "BIN1_CHROMOSOME",
        "BIN2_CHROMOSOME"
      )

    } else {

      candidates <- c(
        "BIN1_CHR",
        "BIN2_CHR",
        "BIN1_CHROMOSOME",
        "BIN2_CHROMOSOME"
      )

      hic_cols <- intersect(
        candidates,
        colnames(hic)
      )
    }

    if (length(hic_cols) == 0) {

      warning(
        "No chromosome columns detected in ",
        basename(f),
        " — skipped."
      )

      output_files[i] <- NA_character_

      next
    }

    # --------------------------------------------------------
    # Remove chrX / chrY interactions
    # --------------------------------------------------------
    xy_regex <- "^(chr)?[XY]$"

    keep <- apply(
      hic[, hic_cols, drop = FALSE],
      1,
      function(x) {

        x <- trimws(
          as.character(x)
        )

        !any(
          grepl(
            xy_regex,
            x,
            ignore.case = TRUE
          ) &
            !is.na(x)
        )
      }
    )

    n_removed <- sum(!keep)

    if (n_removed > 0) {

      message(
        "Removed ",
        n_removed,
        " rows containing chrX/chrY"
      )

    } else {

      message(
        "No chrX/chrY rows found"
      )
    }

    hic_clean <- hic[
      keep,
      ,
      drop = FALSE
    ]

    # --------------------------------------------------------
    # Output cleaned file as CSV
    # --------------------------------------------------------
    out_file <- file.path(
      dirname(f),
      paste0(
        tools::file_path_sans_ext(
          basename(f)
        ),
        "_no_chrX_chrY.csv"
      )
    )

    write.csv(
      hic_clean,
      out_file,
      quote = FALSE,
      row.names = FALSE
    )

    message(
      "Saved: ",
      basename(out_file)
    )

    output_files[i] <- out_file
  }

  invisible(
    output_files[
      !is.na(output_files) &
        output_files != ""
    ]
  )
}
