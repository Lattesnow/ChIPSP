#' Annotate Conventional ChIP and ChIP-SP Regions
#'
#' @description
#' Annotates conventional ChIP-seq peaks and ChIP-SP spatially linked
#' regions to nearby genes using ChIPpeakAnno. Human and mouse genomes
#' are supported.
#'
#' Human reference genomes:
#' \code{hg19}, \code{hg38}
#'
#' Mouse reference genomes:
#' \code{mm9}, \code{mm10}
#'
#' Peaks are annotated using a configurable binding region around genes.
#' Gene symbols are added using the corresponding organism annotation
#' database.
#'
#' Conventional ChIP-seq peaks are ranked by pileup.
#' ChIP-SP regions retain the ChIP-SP score and are ranked by pileup.
#'
#' @param chip_file Character scalar. Path to the conventional ChIP-seq
#'   peak file.
#'
#' @param chipsp_file Character scalar. Path to the ChIP-SP result file.
#'
#' @param species Character scalar specifying species. Must be
#'   \code{"human"} or \code{"mouse"}.
#'
#' @param ref_genome Character scalar specifying reference genome.
#'   Human genomes: \code{"hg19"} or \code{"hg38"}.
#'   Mouse genomes: \code{"mm9"} or \code{"mm10"}.
#'
#' @param binding_bp Numeric scalar specifying the upstream and downstream
#'   annotation distance in base pairs. Default is \code{5000}.
#'
#' @param chip_output Character scalar specifying the output CSV filename
#'   for conventional ChIP annotation.
#'
#' @param chipsp_output Character scalar specifying the output CSV filename
#'   for ChIP-SP annotation.
#'
#' @return Invisibly returns a list containing:
#'
#' \itemize{
#'   \item \code{ChIP}: annotated conventional ChIP-seq table
#'   \item \code{ChIPSP}: annotated ChIP-SP table
#' }
#'
#' @examples
#' \dontrun{
#' annotation_results <- ChIPSPannotation(
#'   chip_file = "Combined_ChIP.bed",
#'   chipsp_file = "ChIPSP_results.csv",
#'   species = "human",
#'   ref_genome = "hg38"
#' )
#' }
#'
#' @export
ChIPSPannotation <- function(
    chip_file,
    chipsp_file,
    species = "human",
    ref_genome = "hg38",
    binding_bp = 5000,
    chip_output = "ChIP_anno_genes_upAdown_UCSC_Control.csv",
    chipsp_output = "ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv"
) {

  # ==========================================================
  # 1. Validate species and genome
  # ==========================================================

  species <- tolower(species)
  ref_genome <- tolower(ref_genome)

  if (!species %in% c("human", "mouse")) {
    stop(
      "`species` must be either \"human\" or \"mouse\"."
    )
  }

    if (species == "human") {
    
      if (!ref_genome %in% c("hg19", "hg38")) {
        stop(
          "For human, `ref_genome` must be \"hg19\" or \"hg38\"."
        )
      }
    
    } else {
    
      if (!ref_genome %in% c("mm9", "mm10")) {
        stop(
          "For mouse, `ref_genome` must be \"mm9\" or \"mm10\"."
        )
      }
    
    }
    
    # ==========================================================
    # 2. Validate files
    # ==========================================================

  if (
    !is.character(chip_file) ||
      length(chip_file) != 1L ||
      is.na(chip_file) ||
      chip_file == ""
  ) {
    stop(
      "`chip_file` must be one valid file path."
    )
  }

  if (!file.exists(chip_file)) {
    stop(
      "Conventional ChIP file not found: ",
      chip_file
    )
  }


  if (
    !is.character(chipsp_file) ||
      length(chipsp_file) != 1L ||
      is.na(chipsp_file) ||
      chipsp_file == ""
  ) {
    stop(
      "`chipsp_file` must be one valid file path."
    )
  }

  if (!file.exists(chipsp_file)) {
    stop(
      "ChIP-SP file not found: ",
      chipsp_file
    )
  }


  # ==========================================================
  # 3. Validate binding distance
  # ==========================================================

  if (
    !is.numeric(binding_bp) ||
      length(binding_bp) != 1L ||
      is.na(binding_bp) ||
      !is.finite(binding_bp) ||
      binding_bp < 0
  ) {
    stop(
      "`binding_bp` must be one finite non-negative numeric value."
    )
  }


  # ==========================================================
  # 4. Check required packages
  # ==========================================================

  required_packages <- c(
    "ChIPpeakAnno",
    "GenomicRanges",
    "IRanges",
    "AnnotationDbi",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "TxDb.Mmusculus.UCSC.mm9.knownGene",
    "TxDb.Mmusculus.UCSC.mm10.knownGene"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "Required package(s) not installed: ",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  }


  # ==========================================================
  # 5. Select genome annotation database
  # ==========================================================

# ==========================================================
# 5. Select genome annotation database
# ==========================================================

if (species == "human") {

  if (ref_genome == "hg19") {

    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

  } else {

    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  }

  orgdb <- "org.Hs.eg.db"

} else {

  if (ref_genome == "mm9") {

    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene

  } else {

    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  }

  orgdb <- "org.Mm.eg.db"
}


  message(
    "Species: ",
    species
  )

  message(
    "Reference genome: ",
    ref_genome
  )

  message(
    "Binding region: ±",
    binding_bp,
    " bp"
  )


  # ==========================================================
  # 6. Build gene annotation GRanges
  # ==========================================================

  annoData <- ChIPpeakAnno::toGRanges(
    txdb,
    feature = "gene"
  )


  # ==========================================================
  # 7. Internal peak-file reader
  # ==========================================================

  read_peak_file <- function(file) {

    extension <- tolower(
      tools::file_ext(file)
    )

    message(
      "Reading: ",
      basename(file)
    )


    # --------------------------------------------------------
    # CSV
    # --------------------------------------------------------

    if (extension == "csv") {

      df <- utils::read.csv(
        file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )


    # --------------------------------------------------------
    # TSV / TXT / TAB
    # --------------------------------------------------------

    } else if (
      extension %in%
        c("tsv", "txt", "tab")
    ) {

      df <- utils::read.delim(
        file,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )


    # --------------------------------------------------------
    # BED
    # --------------------------------------------------------

    } else if (extension == "bed") {

      df <- utils::read.delim(
        file,
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      if (ncol(df) < 3L) {
        stop(
          "BED file must contain at least three columns: ",
          basename(file)
        )
      }

      colnames(df)[1:3] <- c(
        "chr",
        "start",
        "end"
      )

      # Preserve the current ChIP-SP convention:
      # fourth BED column is treated as pileup.
      if (ncol(df) >= 4L) {
        colnames(df)[4] <- "pileup"
      }


    # --------------------------------------------------------
    # XLSX
    # --------------------------------------------------------

    } else if (extension == "xlsx") {

      if (!requireNamespace(
        "readxl",
        quietly = TRUE
      )) {
        stop(
          "Package `readxl` is required to read XLSX files."
        )
      }

      df <- as.data.frame(
        readxl::read_excel(file),
        check.names = FALSE
      )


    # --------------------------------------------------------
    # XLS
    #
    # First try true Excel format.
    # If this fails, try tab-delimited text because older
    # ChIP workflows may use .xls as a text-file extension.
    # --------------------------------------------------------

    } else if (extension == "xls") {

      if (!requireNamespace(
        "readxl",
        quietly = TRUE
      )) {
        stop(
          "Package `readxl` is required to read XLS files."
        )
      }

      df <- tryCatch(

        as.data.frame(
          readxl::read_excel(file),
          check.names = FALSE
        ),

        error = function(e) {

          message(
            "Could not read ",
            basename(file),
            " as native Excel; trying tab-delimited text."
          )

          utils::read.delim(
            file,
            sep = "\t",
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        }
      )


    } else {

      stop(
        "Unsupported peak file type: ",
        extension,
        ". Supported formats are BED, CSV, TSV, TXT, TAB, XLS, and XLSX."
      )
    }


    # --------------------------------------------------------
    # Normalize common column names
    # --------------------------------------------------------

    original_names <- trimws(
      colnames(df)
    )

    upper_names <- toupper(
      original_names
    )

    name_map <- c(
      CHR = "chr",
      CHROM = "chr",
      CHROMOSOME = "chr",
      SEQNAMES = "chr",

      START = "start",
      END = "end",

      PILEUP = "pileup",
      SCORE = "score"
    )

    matched <- upper_names %in%
      names(name_map)

    original_names[
      matched
    ] <- unname(
      name_map[
        upper_names[
          matched
        ]
      ]
    )

    colnames(df) <- original_names


    message(
      "Loaded ",
      nrow(df),
      " rows x ",
      ncol(df),
      " columns."
    )

    return(df)
  }


  # ==========================================================
  # 8. Internal annotation function
  # ==========================================================

  annotate_peak_table <- function(
      df,
      has_score = FALSE,
      out_csv
  ) {

    # --------------------------------------------------------
    # Required columns
    # --------------------------------------------------------

    req <- c(
      "chr",
      "start",
      "end",
      "pileup"
    )

    if (has_score) {
      req <- c(
        req,
        "score"
      )
    }

    missing_columns <- setdiff(
      req,
      colnames(df)
    )

    if (length(missing_columns) > 0L) {
      stop(
        "Required columns missing: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    }


    # --------------------------------------------------------
    # Convert coordinate and signal columns
    # --------------------------------------------------------

    df$chr <- trimws(
      as.character(df$chr)
    )

    df$start <- suppressWarnings(
      as.integer(df$start)
    )

    df$end <- suppressWarnings(
      as.integer(df$end)
    )

    df$pileup <- suppressWarnings(
      as.numeric(df$pileup)
    )

    if (has_score) {
      df$score <- suppressWarnings(
        as.numeric(df$score)
      )
    }


    # --------------------------------------------------------
    # Remove invalid coordinates
    # --------------------------------------------------------

    valid <- (
      !is.na(df$chr) &
        df$chr != "" &
        !is.na(df$start) &
        !is.na(df$end) &
        df$end > df$start &
        !is.na(df$pileup)
    )

    n_removed <- sum(
      !valid
    )

    if (n_removed > 0L) {
      warning(
        "Removed ",
        n_removed,
        " invalid peak rows before annotation."
      )
    }

    df <- df[
      valid,
      ,
      drop = FALSE
    ]


    # --------------------------------------------------------
    # Generate GRanges
    # --------------------------------------------------------

    gr <- GenomicRanges::GRanges(
      seqnames = df$chr,
      ranges = IRanges::IRanges(
        start = df$start,
        end = df$end
      ),
      pileup = df$pileup
    )

    if (has_score) {

      S4Vectors::mcols(
        gr
      )$score <- df$score
    }


    # --------------------------------------------------------
    # Annotate peaks
    # --------------------------------------------------------

    anno <- ChIPpeakAnno::annotatePeakInBatch(
      gr,
      AnnotationData = annoData,
      output = "both",
      bindingRegion = c(
        -binding_bp,
        binding_bp
      )
    )


    # --------------------------------------------------------
    # Add gene symbols
    # --------------------------------------------------------

    anno <- ChIPpeakAnno::addGeneIDs(
      anno,
      orgdb,
      IDs2Add = "symbol",
      feature_id_type = "entrez_id"
    )


    anno_df <- as.data.frame(
      anno
    )


    # --------------------------------------------------------
    # Merge original pileup and score values
    # --------------------------------------------------------

    gr_df <- as.data.frame(
      gr
    )

    keep_columns <- c(
      "seqnames",
      "start",
      "end",
      "pileup"
    )

    if (has_score) {
      keep_columns <- c(
        keep_columns,
        "score"
      )
    }

    gr_df <- gr_df[
      ,
      keep_columns,
      drop = FALSE
    ]


    anno_df <- merge(
      anno_df,
      gr_df,
      by = c(
        "seqnames",
        "start",
        "end"
      ),
      all.x = TRUE
    )


    # --------------------------------------------------------
    # Detect pileup column after merge
    # --------------------------------------------------------

    pileup_candidates <- c(
      "pileup",
      "pileup.x",
      "pileup.y"
    )

    pileup_col <- intersect(
      pileup_candidates,
      colnames(anno_df)
    )

    if (length(pileup_col) == 0L) {
      stop(
        "No pileup column found after annotation merge."
      )
    }

    pileup_col <- pileup_col[1]


    # --------------------------------------------------------
    # Detect score column for ChIP-SP
    # --------------------------------------------------------

    if (has_score) {

      score_candidates <- c(
        "score",
        "score.x",
        "score.y"
      )

      score_col <- intersect(
        score_candidates,
        colnames(anno_df)
      )

      if (length(score_col) > 0L) {

        score_col <- score_col[1]

        if (score_col != "score") {
            anno_df$score <- anno_df[[score_col]]
        }
      }
    }


    # --------------------------------------------------------
    # Remove entries without gene symbols
    # --------------------------------------------------------

    if (!"symbol" %in% colnames(anno_df)) {
      stop(
        "Gene symbol column was not generated during annotation."
      )
    }

    anno_df <- anno_df[
      !is.na(anno_df$symbol) &
        anno_df$symbol != "",
      ,
      drop = FALSE
    ]


    # --------------------------------------------------------
    # Sort by pileup
    # --------------------------------------------------------

    anno_df <- anno_df[
      order(
        -anno_df[[pileup_col]]
      ),
      ,
      drop = FALSE
    ]


    # --------------------------------------------------------
    # Dense rank by pileup
    # --------------------------------------------------------

    anno_df$rank <- match(
        anno_df[[pileup_col]],
      sort(
        unique(
         anno_df[[pileup_col]]
        ),
        decreasing = TRUE
      )
    )


    # --------------------------------------------------------
    # Write annotation table
    # --------------------------------------------------------

    utils::write.csv(
      anno_df,
      file = out_csv,
      row.names = FALSE,
      quote = FALSE
    )

    message(
      "Wrote: ",
      normalizePath(
        out_csv,
        mustWork = FALSE
      )
    )

    invisible(
      anno_df
    )
  }


  # ==========================================================
  # 9. Read conventional ChIP peaks
  # ==========================================================

  message(
    "--------------------------------------------------"
  )

  message(
    "Annotating conventional ChIP-seq peaks"
  )

  message(
    "--------------------------------------------------"
  )

  chip <- read_peak_file(
    chip_file
  )


  chip_annotation <- annotate_peak_table(
    df = chip,
    has_score = FALSE,
    out_csv = chip_output
  )


  # ==========================================================
  # 10. Read ChIP-SP output
  # ==========================================================

  message(
    "--------------------------------------------------"
  )

  message(
    "Annotating ChIP-SP regions"
  )

  message(
    "--------------------------------------------------"
  )

  chipsp <- read_peak_file(
    chipsp_file
  )


  chipsp_annotation <- annotate_peak_table(
    df = chipsp,
    has_score = TRUE,
    out_csv = chipsp_output
  )


  # ==========================================================
  # 11. Summary
  # ==========================================================

  message(
    "--------------------------------------------------"
  )

  message(
    "ChIP-SP annotation completed."
  )

  message(
    "Conventional ChIP annotated rows: ",
    nrow(
      chip_annotation
    )
  )

  message(
    "ChIP-SP annotated rows: ",
    nrow(
      chipsp_annotation
    )
  )

  message(
    "--------------------------------------------------"
  )


  invisible(
    list(
      ChIP = chip_annotation,
      ChIPSP = chipsp_annotation
    )
  )
}
