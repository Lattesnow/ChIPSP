#' Run Enrichr Analysis for ChIP and ChIP-SP Gene Sets
#'
#' @description
#' Performs Enrichr pathway enrichment analysis for genes annotated from
#' conventional ChIP-seq and ChIP-SP results.
#'
#' Five gene sets are generated:
#' \itemize{
#'   \item ChIP
#'   \item ChIPSP
#'   \item ChIPSP_only
#'   \item ChIP_only
#'   \item Shared
#' }
#'
#' Enrichr result tables are written for every selected database.
#' Dotplots are generated for ChIP, ChIPSP, ChIPSP_only, and Shared genes.
#'
#' Pathways are ranked by nominal P value.
#' Dot color represents -log10(P value).
#' Dot size represents the number of overlapping genes.
#' The x-axis represents the gene ratio:
#' overlapping genes / total genes in the pathway.
#'
#' @param chip_file Character scalar. CSV file containing annotated
#'   conventional ChIP-seq genes. Default is
#'   \code{"ChIP_anno_genes_upAdown_UCSC_Control.csv"}.
#'
#' @param chipsp_file Character scalar. CSV file containing annotated
#'   ChIP-SP genes. Default is
#'   \code{"ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv"}.
#'
#' @param out_dir Character scalar. Output directory.
#'   Default is \code{"Enrichr_dotplots"}.
#'
#' @param databases Character vector containing Enrichr database names.
#'
#' @param symbol_col Character scalar specifying the gene-symbol column.
#'   Default is \code{"symbol"}.
#'
#' @param top_n Integer specifying the number of pathways shown in each
#'   dotplot. Default is \code{10}.
#'
#' @param enrichment_top_n Integer specifying the number of pathways retained
#'   from each Enrichr database when generating the combined result table.
#'   Default is \code{50}.
#'
#' @param p_cutoff Numeric scalar specifying the nominal P-value cutoff for
#'   dotplots. Default is \code{0.05}.
#'
#' @param uppercase_genes Logical. If TRUE, gene symbols are converted to
#'   uppercase before Enrichr analysis. Default is TRUE.
#'
#' @return Invisibly returns a list containing gene sets, gene summary,
#'   Enrichr results, combined enrichment results, and generated plots.
#'
#' @examples
#' \dontrun{
#' enrichment_results <- ChIPSPenrichment()
#' }
#'
#' @export
ChIPSPenrichment <- function(
    chip_file = "ChIP_anno_genes_upAdown_UCSC_Control.csv",
    chipsp_file = "ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv",
    out_dir = "Enrichr_dotplots",
    databases = c(
      "GO_Biological_Process_2023",
      "GO_Cellular_Component_2023",
      "GO_Molecular_Function_2023",
      "KEGG_2021_Human",
      "Reactome_2022",
      "ChEA_2022",
      "MSigDB_Hallmark_2020"
    ),
    symbol_col = "symbol",
    top_n = 10,
    enrichment_top_n = 50,
    p_cutoff = 0.05,
    uppercase_genes = TRUE
) {

  # ==========================================================
  # 1. Check required packages
  # ==========================================================

  required_packages <- c(
    "enrichR",
    "ggplot2",
    "stringr",
    "scales"
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
  # 2. Validate input files
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
      "ChIP annotation file not found: ",
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
      "ChIP-SP annotation file not found: ",
      chipsp_file
    )
  }


  # ==========================================================
  # 3. Validate arguments
  # ==========================================================

  if (
    !is.character(out_dir) ||
      length(out_dir) != 1L ||
      is.na(out_dir) ||
      out_dir == ""
  ) {
    stop(
      "`out_dir` must be one valid directory path."
    )
  }


  if (
    !is.character(databases) ||
      length(databases) == 0L
  ) {
    stop(
      "`databases` must contain at least one Enrichr database."
    )
  }


  if (
    !is.character(symbol_col) ||
      length(symbol_col) != 1L ||
      is.na(symbol_col) ||
      symbol_col == ""
  ) {
    stop(
      "`symbol_col` must be one valid column name."
    )
  }


  if (
    !is.numeric(top_n) ||
      length(top_n) != 1L ||
      is.na(top_n) ||
      top_n < 1
  ) {
    stop(
      "`top_n` must be a positive integer."
    )
  }

  top_n <- as.integer(top_n)


  if (
    !is.numeric(enrichment_top_n) ||
      length(enrichment_top_n) != 1L ||
      is.na(enrichment_top_n) ||
      enrichment_top_n < 1
  ) {
    stop(
      "`enrichment_top_n` must be a positive integer."
    )
  }

  enrichment_top_n <- as.integer(
    enrichment_top_n
  )


  if (
    !is.numeric(p_cutoff) ||
      length(p_cutoff) != 1L ||
      is.na(p_cutoff) ||
      !is.finite(p_cutoff) ||
      p_cutoff <= 0 ||
      p_cutoff > 1
  ) {
    stop(
      "`p_cutoff` must be greater than 0 and less than or equal to 1."
    )
  }


  if (
    !is.logical(uppercase_genes) ||
      length(uppercase_genes) != 1L ||
      is.na(uppercase_genes)
  ) {
    stop(
      "`uppercase_genes` must be TRUE or FALSE."
    )
  }


  # ==========================================================
  # 4. Create output directory
  # ==========================================================

  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )


  # ==========================================================
  # 5. Initialize Enrichr
  #
  # Important:
  # enrichR::setEnrichrSite() can fail when enrichR is used
  # through :: without first attaching the package because
  # enrichR.sites.base.address may not have been initialized.
  #
  # Set the required Enrichr options directly.
  # ==========================================================

  options(
    enrichR.sites.base.address =
      "https://maayanlab.cloud/"
  )

  options(
    enrichR.base.address =
      "https://maayanlab.cloud/Enrichr/"
  )

  options(
    speedrichr.base.address =
      "https://maayanlab.cloud/speedrichr/api/"
  )

  options(
    enrichR.live = TRUE
  )

  options(
    enrichR.quiet = FALSE
  )

  options(
    enrichR.sites = c(
      "Enrichr",
      "FlyEnrichr",
      "WormEnrichr",
      "YeastEnrichr",
      "FishEnrichr",
      "OxEnrichr"
    )
  )


  message(
    "Enrichr server: ",
    getOption("enrichR.base.address")
  )


  # ==========================================================
  # 6. Reproducibility
  # ==========================================================

  set.seed(
    1234
  )


  # ==========================================================
  # 7. Internal function:
  # Read gene list
  # ==========================================================

  get_gene_list <- function(
      file,
      symbol_column
  ) {

    tab <- utils::read.csv(
      file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )


    if (!symbol_column %in% colnames(tab)) {
      stop(
        "Column `",
        symbol_column,
        "` was not found in ",
        basename(file),
        "."
      )
    }


    genes <- unique(
      as.character(
        tab[[symbol_column]]
      )
    )


    genes <- trimws(
      genes
    )


    genes <- genes[
      !is.na(genes) &
        genes != ""
    ]


    if (uppercase_genes) {
      genes <- toupper(
        genes
      )
    }


    genes <- unique(
      genes
    )


    return(
      genes
    )
  }


  # ==========================================================
  # 8. Internal function:
  # Run Enrichr and save result tables
  # ==========================================================

  run_enrichr_save <- function(
      gene_list,
      group_name
  ) {

    message(
      "--------------------------------------------------"
    )

    message(
      "Running Enrichr for: ",
      group_name
    )

    message(
      "Input genes: ",
      length(gene_list)
    )


    group_dir <- file.path(
      out_dir,
      group_name
    )


    dir.create(
      group_dir,
      showWarnings = FALSE,
      recursive = TRUE
    )


    if (length(gene_list) < 5L) {

      warning(
        group_name,
        " has fewer than 5 genes. Skipping Enrichr analysis."
      )

      return(
        NULL
      )
    }


    enr <- enrichR::enrichr(
      gene_list,
      databases
    )


    if (
      is.null(enr) ||
        length(enr) == 0L
    ) {

      warning(
        "No Enrichr results returned for ",
        group_name,
        "."
      )

      return(
        NULL
      )
    }


    for (db in names(enr)) {

      result_file <- file.path(
        group_dir,
        paste0(
          db,
          "_",
          group_name,
          ".csv"
        )
      )


      utils::write.csv(
        enr[[db]],
        file = result_file,
        row.names = FALSE
      )
    }


    return(
      enr
    )
  }


  # ==========================================================
  # 9. Internal function:
  # Extract top Enrichr pathways
  # ==========================================================

  extract_top_terms <- function(
      enr_list,
      group_name,
      keep_n
  ) {

    if (is.null(enr_list)) {
      return(
        NULL
      )
    }


    result_list <- list()


    for (db in names(enr_list)) {

      tab <- as.data.frame(
        enr_list[[db]],
        stringsAsFactors = FALSE
      )


      if (nrow(tab) == 0L) {
        next
      }


      required_columns <- c(
        "Term",
        "P.value",
        "Adjusted.P.value",
        "Overlap"
      )


      missing_columns <- setdiff(
        required_columns,
        colnames(tab)
      )


      if (length(missing_columns) > 0L) {

        warning(
          "Skipping database ",
          db,
          " because required column(s) are missing: ",
          paste(
            missing_columns,
            collapse = ", "
          )
        )

        next
      }


      tab$Database <- db
      tab$Group <- group_name
      tab$Term_clean <- as.character(
        tab$Term
      )


      tab$P.value <- suppressWarnings(
        as.numeric(
          tab$P.value
        )
      )


      tab$Adjusted.P.value <- suppressWarnings(
        as.numeric(
          tab$Adjusted.P.value
        )
      )


      tab$neg_log10_p <- -log10(
        tab$P.value
      )


      tab$neg_log10_adjP <- -log10(
        tab$Adjusted.P.value
      )


      tab$overlap_hit <- suppressWarnings(
        as.numeric(
          sub(
            "/.*",
            "",
            tab$Overlap
          )
        )
      )


      tab$overlap_total <- suppressWarnings(
        as.numeric(
          sub(
            ".*/",
            "",
            tab$Overlap
          )
        )
      )


      tab$GeneRatio <- (
        tab$overlap_hit /
          tab$overlap_total
      )


      tab <- tab[
        !is.na(tab$P.value),
        ,
        drop = FALSE
      ]


      if (nrow(tab) == 0L) {
        next
      }


      tab <- tab[
        order(
          tab$P.value,
          decreasing = FALSE,
          na.last = TRUE
        ),
        ,
        drop = FALSE
      ]


      if (nrow(tab) > keep_n) {

        tab <- tab[
          seq_len(
            keep_n
          ),
          ,
          drop = FALSE
        ]
      }


      result_list[[db]] <- tab
    }


    if (length(result_list) == 0L) {
      return(
        NULL
      )
    }


    out <- do.call(
      rbind,
      result_list
    )


    rownames(out) <- NULL


    return(
      out
    )
  }


  # ==========================================================
  # 10. Internal function:
  # Nature-style theme
  # ==========================================================

  theme_nature_dot <- function(
      base_size = 8
  ) {

    ggplot2::theme_classic(
      base_size = base_size,
      base_family = "Arial"
    ) +

      ggplot2::theme(

        text = ggplot2::element_text(
          family = "Arial",
          color = "black"
        ),

        plot.title = ggplot2::element_text(
          family = "Arial",
          size = 10,
          face = "bold",
          hjust = 0
        ),

        axis.title = ggplot2::element_text(
          family = "Arial",
          size = 9
        ),

        axis.text = ggplot2::element_text(
          family = "Arial",
          size = 8,
          color = "black"
        ),

        axis.text.x = ggplot2::element_text(
          family = "Arial",
          size = 8,
          angle = 45,
          hjust = 1
        ),

        axis.text.y = ggplot2::element_text(
          family = "Arial",
          size = 10,
          color = "black"
        ),

        axis.line = ggplot2::element_line(
          linewidth = 0.5,
          color = "black"
        ),

        axis.ticks = ggplot2::element_line(
          linewidth = 0.5,
          color = "black"
        ),

        panel.grid.major = ggplot2::element_line(
          linewidth = 0.25,
          color = "grey88"
        ),

        panel.grid.minor = ggplot2::element_blank(),

        legend.title = ggplot2::element_text(
          family = "Arial",
          size = 8
        ),

        legend.text = ggplot2::element_text(
          family = "Arial",
          size = 7
        ),

        plot.margin = ggplot2::margin(
          5,
          8,
          5,
          5
        )
      )
  }


  # ==========================================================
  # 11. Internal function:
  # Enrichr dotplot
  #
  # Ranking:
  # P value, smallest first
  #
  # Dot color:
  # -log10(P value)
  #
  # Dot size:
  # overlapping gene count
  #
  # X-axis:
  # overlap_hit / overlap_total
  # ==========================================================

  plot_enrichr_dotplot <- function(
      summary_table,
      group_name,
      database_name
  ) {

    if (
      is.null(summary_table) ||
        nrow(summary_table) == 0L
    ) {

      return(
        NULL
      )
    }


    keep <- (
      summary_table$Group == group_name &
        summary_table$Database == database_name &
        !is.na(summary_table$P.value) &
        summary_table$P.value < p_cutoff
    )


    df <- summary_table[
      keep,
      ,
      drop = FALSE
    ]


    if (nrow(df) == 0L) {

      message(
        "No nominal P < ",
        p_cutoff,
        " terms for ",
        group_name,
        " - ",
        database_name
      )

      return(
        NULL
      )
    }


    # --------------------------------------------------------
    # Rank pathways by P value
    # --------------------------------------------------------

    df <- df[
      order(
        df$P.value,
        decreasing = FALSE
      ),
      ,
      drop = FALSE
    ]


    if (nrow(df) > top_n) {

      df <- df[
        seq_len(
          top_n
        ),
        ,
        drop = FALSE
      ]
    }


    # --------------------------------------------------------
    # Clean pathway names
    # --------------------------------------------------------

    df$Term_clean <- stringr::str_replace(
      df$Term_clean,
      "Homo sapiens",
      ""
    )


    df$Term_clean <- stringr::str_replace_all(
      df$Term_clean,
      "_",
      " "
    )


    df$Term_clean <- stringr::str_replace_all(
      df$Term_clean,
      "\\s+",
      " "
    )


    df$Term_clean <- stringr::str_replace(
      df$Term_clean,
      "\\s*R-HSA-[0-9]+",
      ""
    )


    df$Term_clean <- stringr::str_trim(
      df$Term_clean
    )


    df$Term_clean <- stringr::str_trunc(
      df$Term_clean,
      58
    )


    # --------------------------------------------------------
    # Smallest P value at top
    # --------------------------------------------------------

    df$Term_clean <- factor(
      df$Term_clean,
      levels = rev(
        unique(
          df$Term_clean
        )
      )
    )


    # --------------------------------------------------------
    # Dot color
    # --------------------------------------------------------

    df$log10_p <- -log10(
      df$P.value
    )


    # --------------------------------------------------------
    # Remove rows with invalid plot values
    # --------------------------------------------------------

    valid_plot <- (
      !is.na(df$GeneRatio) &
        is.finite(df$GeneRatio) &
        !is.na(df$overlap_hit) &
        is.finite(df$overlap_hit) &
        !is.na(df$log10_p) &
        is.finite(df$log10_p)
    )


    df <- df[
      valid_plot,
      ,
      drop = FALSE
    ]


    if (nrow(df) == 0L) {

      message(
        "No valid plotting values for ",
        group_name,
        " - ",
        database_name
      )

      return(
        NULL
      )
    }


    max_ratio <- max(
      df$GeneRatio,
      na.rm = TRUE
    )


    max_log10_p <- max(
      df$log10_p,
      na.rm = TRUE
    )


    # Since plotted pathways have P < 0.05,
    # max_log10_p should normally be > 1.
    # This protects against an invalid color-scale range.
    color_upper <- max(
      1.01,
      max_log10_p
    )


    if (
      !is.finite(max_ratio) ||
        max_ratio <= 0
    ) {

      max_ratio <- 1
    }


    # ========================================================
    # Dotplot
    # ========================================================

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = GeneRatio,
        y = Term_clean
      )
    ) +

      ggplot2::geom_point(
        ggplot2::aes(
          size = overlap_hit,
          color = log10_p
        ),
        alpha = 0.95
      ) +

      ggplot2::scale_color_gradientn(
        colors = c(
          "#2F8DD8",
          "#9B8FBF",
          "#E46A6A"
        ),
        limits = c(
          1,
          color_upper
        ),
        oob = scales::squish,
        name = expression(
          -log[10](P)
        )
      ) +

      ggplot2::scale_size_continuous(
        range = c(
          1.8,
          6.5
        ),
        breaks = scales::pretty_breaks(
          n = 4
        ),
        name = "Gene count"
      ) +

      ggplot2::scale_x_continuous(
        limits = c(
          0,
          max_ratio * 1.12
        ),
        expand = ggplot2::expansion(
          mult = c(
            0.03,
            0.08
          )
        )
      ) +

      ggplot2::labs(
        title = paste0(
          group_name,
          " - ",
          database_name
        ),
        x = "Gene ratio",
        y = NULL
      ) +

      theme_nature_dot() +

      ggplot2::theme(
        plot.margin = ggplot2::margin(
          5,
          35,
          5,
          20
        ),
        legend.position = "right"
      )


    # ========================================================
    # Output filenames
    # ========================================================

    safe_db <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      database_name
    )


    safe_group <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      group_name
    )


    pdf_file <- file.path(
      out_dir,
      paste0(
        "dotplot_",
        safe_group,
        "_",
        safe_db,
        ".pdf"
      )
    )


    png_file <- file.path(
      out_dir,
      paste0(
        "dotplot_",
        safe_group,
        "_",
        safe_db,
        ".png"
      )
    )


    # ========================================================
    # Save PDF
    # ========================================================

    if (capabilities("cairo")) {

      ggplot2::ggsave(
        filename = pdf_file,
        plot = p,
        device = grDevices::cairo_pdf,
        width = 6.6,
        height = 3.6,
        units = "in"
      )

    } else {

      warning(
        "Cairo graphics are not available. ",
        "Using the standard PDF device for ",
        basename(pdf_file),
        "."
      )

      ggplot2::ggsave(
        filename = pdf_file,
        plot = p,
        device = "pdf",
        width = 6.6,
        height = 3.6,
        units = "in"
      )
    }


    # ========================================================
    # Save PNG
    # ========================================================

    ggplot2::ggsave(
      filename = png_file,
      plot = p,
      width = 6.6,
      height = 3.6,
      units = "in",
      dpi = 600,
      bg = "white"
    )


    return(
      p
    )
  }


  # ==========================================================
  # 12. Read input gene lists
  # ==========================================================

  message(
    "--------------------------------------------------"
  )

  message(
    "Reading annotated gene lists"
  )

  message(
    "--------------------------------------------------"
  )


  chip_genes <- get_gene_list(
    file = chip_file,
    symbol_column = symbol_col
  )


  chipsp_genes <- get_gene_list(
    file = chipsp_file,
    symbol_column = symbol_col
  )


  chipsp_only_genes <- setdiff(
    chipsp_genes,
    chip_genes
  )


  chip_only_genes <- setdiff(
    chip_genes,
    chipsp_genes
  )


  shared_genes <- intersect(
    chip_genes,
    chipsp_genes
  )


  # ==========================================================
  # 13. Gene-set summary
  # ==========================================================

  gene_summary <- data.frame(

    Category = c(
      "ChIP",
      "ChIPSP",
      "ChIPSP_only",
      "ChIP_only",
      "Shared"
    ),

    Gene_count = c(
      length(chip_genes),
      length(chipsp_genes),
      length(chipsp_only_genes),
      length(chip_only_genes),
      length(shared_genes)
    ),

    stringsAsFactors = FALSE
  )


  utils::write.csv(
    gene_summary,
    file = file.path(
      out_dir,
      "Gene_set_summary.csv"
    ),
    row.names = FALSE
  )


  writeLines(
    chip_genes,
    con = file.path(
      out_dir,
      "ChIP_gene_list.txt"
    )
  )


  writeLines(
    chipsp_genes,
    con = file.path(
      out_dir,
      "ChIPSP_gene_list.txt"
    )
  )


  writeLines(
    chipsp_only_genes,
    con = file.path(
      out_dir,
      "ChIPSP_only_gene_list.txt"
    )
  )


  writeLines(
    chip_only_genes,
    con = file.path(
      out_dir,
      "ChIP_only_gene_list.txt"
    )
  )


  writeLines(
    shared_genes,
    con = file.path(
      out_dir,
      "Shared_gene_list.txt"
    )
  )


  print(
    gene_summary
  )


  # ==========================================================
  # 14. Run Enrichr
  # ==========================================================

  enr_chip <- run_enrichr_save(
    gene_list = chip_genes,
    group_name = "ChIP"
  )


  enr_chipsp <- run_enrichr_save(
    gene_list = chipsp_genes,
    group_name = "ChIPSP"
  )


  enr_chipsp_only <- run_enrichr_save(
    gene_list = chipsp_only_genes,
    group_name = "ChIPSP_only"
  )


  enr_chip_only <- run_enrichr_save(
    gene_list = chip_only_genes,
    group_name = "ChIP_only"
  )


  enr_shared <- run_enrichr_save(
    gene_list = shared_genes,
    group_name = "Shared"
  )


  # ==========================================================
  # 15. Extract and combine enrichment results
  # ==========================================================

  summary_parts <- list(

    extract_top_terms(
      enr_list = enr_chip,
      group_name = "ChIP",
      keep_n = enrichment_top_n
    ),

    extract_top_terms(
      enr_list = enr_chipsp,
      group_name = "ChIPSP",
      keep_n = enrichment_top_n
    ),

    extract_top_terms(
      enr_list = enr_chipsp_only,
      group_name = "ChIPSP_only",
      keep_n = enrichment_top_n
    ),

    extract_top_terms(
      enr_list = enr_chip_only,
      group_name = "ChIP_only",
      keep_n = enrichment_top_n
    ),

    extract_top_terms(
      enr_list = enr_shared,
      group_name = "Shared",
      keep_n = enrichment_top_n
    )
  )


  summary_parts <- Filter(
    Negate(
      is.null
    ),
    summary_parts
  )


  if (length(summary_parts) == 0L) {

    warning(
      "No Enrichr results were available for any gene group."
    )

    summary_all <- data.frame()

  } else {

    summary_all <- do.call(
      rbind,
      summary_parts
    )

    rownames(
      summary_all
    ) <- NULL
  }


  utils::write.csv(
    summary_all,
    file = file.path(
      out_dir,
      "Combined_Enrichr_results_all_groups.csv"
    ),
    row.names = FALSE
  )


  # ==========================================================
  # 16. Generate dotplots
  # ==========================================================

  groups_to_plot <- c(
    "ChIP",
    "ChIPSP",
    "ChIPSP_only",
    "Shared"
  )


  plots <- list()


  if (nrow(summary_all) > 0L) {

    for (g in groups_to_plot) {

      plots[[g]] <- list()


      for (db in databases) {

        p <- plot_enrichr_dotplot(
          summary_table = summary_all,
          group_name = g,
          database_name = db
        )


        plots[[g]][[db]] <- p
      }
    }
  }


  # ==========================================================
  # 17. Completion summary
  # ==========================================================

  message(
    "--------------------------------------------------"
  )

  message(
    "ChIP-SP Enrichr analysis completed."
  )

  message(
    "Output directory: ",
    normalizePath(
      out_dir,
      mustWork = FALSE
    )
  )

  message(
    "ChIP genes: ",
    length(chip_genes)
  )

  message(
    "ChIP-SP genes: ",
    length(chipsp_genes)
  )

  message(
    "ChIP-SP-only genes: ",
    length(chipsp_only_genes)
  )

  message(
    "ChIP-only genes: ",
    length(chip_only_genes)
  )

  message(
    "Shared genes: ",
    length(shared_genes)
  )

  message(
    "--------------------------------------------------"
  )


  # ==========================================================
  # 18. Return
  # ==========================================================

  invisible(
    list(

      gene_summary = gene_summary,

      genes = list(
        ChIP = chip_genes,
        ChIPSP = chipsp_genes,
        ChIPSP_only = chipsp_only_genes,
        ChIP_only = chip_only_genes,
        Shared = shared_genes
      ),

      enrichment = list(
        ChIP = enr_chip,
        ChIPSP = enr_chipsp,
        ChIPSP_only = enr_chipsp_only,
        ChIP_only = enr_chip_only,
        Shared = enr_shared
      ),

      combined_results = summary_all,

      plots = plots
    )
  )
}
