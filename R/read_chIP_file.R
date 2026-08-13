#' Automatically find and read a ChIP-seq peak file
#'
#' Searches for a file whose name ends with "ChIP" followed by a
#' supported extension and reads it into a data.table.
#'
#' Supported formats are BED, CSV, TSV, TXT, TAB, XLS, and XLSX.
#'
#' If `file` is not supplied, the function searches `path`
#' automatically for files matching:
#'
#'   *ChIP.bed
#'   *ChIP.csv
#'   *ChIP.tsv
#'   *ChIP.txt
#'   *ChIP.tab
#'   *ChIP.xls
#'   *ChIP.xlsx
#'
#' @param file Optional path to one ChIP-seq file.
#' @param path Directory searched when `file = NULL`.
#'   Default is the current working directory.
#'
#' @return A data.table containing the ChIP-seq peaks.
#'
#' @export
readChIPFile <- function(
    file = NULL,
    path = getwd()
) {

  # ----------------------------------------------------------
  # 1. Automatically locate ChIP file
  # ----------------------------------------------------------

  if (is.null(file)) {

    chip_files <- list.files(
      path = path,
      pattern = "ChIP\\.(bed|csv|tsv|txt|tab|xls|xlsx)$",
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(chip_files) == 0L) {
      stop(
        "No ChIP file found in: ",
        normalizePath(path, mustWork = FALSE),
        "\nExpected a filename ending with ",
        "ChIP.bed, ChIP.csv, ChIP.tsv, ChIP.txt, ",
        "ChIP.tab, ChIP.xls, or ChIP.xlsx."
      )
    }

    if (length(chip_files) > 1L) {

      stop(
        "More than one ChIP file was found:\n",
        paste(
          basename(chip_files),
          collapse = "\n"
        ),
        "\n\nPlease specify the desired file using `file = ...`."
      )
    }

    file <- chip_files[1]
  }


  # ----------------------------------------------------------
  # 2. Validate supplied file
  # ----------------------------------------------------------

  if (
    !is.character(file) ||
      length(file) != 1L ||
      is.na(file) ||
      file == ""
  ) {
    stop("`file` must contain one valid ChIP file path.")
  }

  if (!file.exists(file)) {
    stop(
      "ChIP file not found: ",
      file
    )
  }


  message(
    "ChIP file detected: ",
    basename(file)
  )


  # ----------------------------------------------------------
  # 3. Determine file extension
  # ----------------------------------------------------------

  extension <- tolower(
    tools::file_ext(file)
  )


  # ----------------------------------------------------------
  # 4. Read file according to format
  # ----------------------------------------------------------

  if (extension == "csv") {

    chip <- data.table::fread(
      file,
      sep = ",",
      header = TRUE,
      fill = TRUE,
      data.table = TRUE
    )


  } else if (extension %in% c(
    "tsv",
    "txt",
    "tab"
  )) {

    chip <- data.table::fread(
      file,
      sep = "\t",
      header = TRUE,
      fill = TRUE,
      data.table = TRUE
    )


  } else if (extension == "bed") {

    chip <- data.table::fread(
      file,
      sep = "\t",
      header = FALSE,
      fill = TRUE,
      data.table = TRUE
    )

    if (ncol(chip) < 3L) {
      stop(
        "BED file must contain at least three columns."
      )
    }

    # BED:
    # column 1 = chromosome
    # column 2 = start
    # column 3 = end
    data.table::setnames(
      chip,
      old = names(chip)[1:3],
      new = c(
        "chr",
        "start",
        "end"
      )
    )

    # Preserve current ChIP-SP behavior:
    # if a 4th column exists, use it as pileup
    if (ncol(chip) >= 4L) {

      data.table::setnames(
        chip,
        old = names(chip)[4],
        new = "pileup"
      )
    }


  } else if (extension == "xlsx") {

    if (!requireNamespace(
      "readxl",
      quietly = TRUE
    )) {
      stop(
        "Package `readxl` is required ",
        "to read Excel files."
      )
    }

    chip <- data.table::as.data.table(
      readxl::read_excel(file)
    )


  } else if (extension == "xls") {

    # --------------------------------------------------------
    # XLS is handled specially.
    #
    # First try as a real Excel XLS file.
    #
    # If that fails, try as tab-delimited text. This supports
    # BED/text files that were renamed with an .xls extension.
    # --------------------------------------------------------

    if (!requireNamespace(
      "readxl",
      quietly = TRUE
    )) {
      stop(
        "Package `readxl` is required ",
        "to read Excel files."
      )
    }

    chip <- tryCatch(

      {

        message(
          "Attempting to read .xls as an Excel file..."
        )

        data.table::as.data.table(
          readxl::read_excel(file)
        )
      },

      error = function(e) {

        message(
          "Excel reading failed. ",
          "Trying tab-delimited text format..."
        )

        data.table::fread(
          file,
          sep = "\t",
          header = FALSE,
          fill = TRUE,
          data.table = TRUE
        )
      }
    )


    # --------------------------------------------------------
    # If fallback produced BED-like unnamed columns,
    # assign chr/start/end
    # --------------------------------------------------------

    generic_names <- all(
      grepl(
        "^V[0-9]+$",
        names(chip)
      )
    )

    if (
      generic_names &&
        ncol(chip) >= 3L
    ) {

      data.table::setnames(
        chip,
        old = names(chip)[1:3],
        new = c(
          "chr",
          "start",
          "end"
        )
      )

      if (ncol(chip) >= 4L) {

        data.table::setnames(
          chip,
          old = names(chip)[4],
          new = "pileup"
        )
      }
    }


  } else {

    stop(
      "Unsupported ChIP file type: ",
      extension,
      ". Supported formats are ",
      "bed, csv, tsv, txt, tab, xls, and xlsx."
    )
  }


  # ----------------------------------------------------------
  # 5. Report result
  # ----------------------------------------------------------

  message(
    "Loaded ChIP file: ",
    basename(file),
    " [",
    nrow(chip),
    " rows x ",
    ncol(chip),
    " columns]"
  )


  attr(
    chip,
    "chip_file"
  ) <- file


  return(chip)
}
