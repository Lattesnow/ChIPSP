#' Run Enrichr Analysis for ChIP and ChIP-SP Gene Sets
#'
#' @description
#' Performs Enrichr analysis for genes annotated from conventional ChIP-seq
#' and ChIP-SP results. The function generates ChIP, ChIP-SP, ChIP-SP-only,
#' ChIP-only, and shared gene sets, runs enrichment analysis, saves enrichment
#' tables, and creates manuscript-style dot plots.
#'
#' Pathways are ranked by nominal P value. Dot color represents
#' -log10(P value), dot size represents overlapping gene count, and the
#' x-axis represents gene ratio (overlapping genes / total pathway genes).
#'
#' @param chip_file Character scalar. CSV file containing annotated
#'   conventional ChIP-seq genes. Default is
#'   \code{"ChIP_anno_genes_upAdown_UCSC_Control.csv"}.
#' @param chipsp_file Character scalar. CSV file containing annotated
#'   ChIP-SP genes. Default is
#'   \code{"ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv"}.
#' @param out_dir Character scalar. Output directory. Default is
#'   \code{"Enrichr_dotplots"}.
#' @param symbol_col Character scalar. Gene-symbol column in the input files.
#'   Default is \code{"symbol"}.
#' @param databases Character vector of Enrichr database names.
#' @param top_n_summary Integer. Maximum number of terms retained per database
#'   in the combined result table. Default is 50.
#' @param top_n_plot Integer. Maximum number of pathways shown per dot plot.
#'   Default is 10.
#' @param p_cutoff Numeric scalar. Nominal P-value cutoff used for plotting.
#'   Default is 0.05.
#' @param groups_to_plot Character vector specifying groups for which plots
#'   are generated. Default is ChIP, ChIPSP, ChIPSP_only, and Shared.
#' @param uppercase_genes Logical. If TRUE, gene symbols are converted to
#'   uppercase before Enrichr analysis, matching the original workflow.
#'   Default is TRUE.
#' @param seed Integer random seed. Default is 1234.
#'
#' @return Invisibly returns a list containing gene sets, the gene summary,
#'   individual Enrichr results, the combined enrichment table, and plots.
#'
#' @examples
#' \dontrun{
#' enrich_results <- ChIPSPenrichment()
#' }
#'
#' @export
ChIPSPenrichment <- function(
    chip_file = "ChIP_anno_genes_upAdown_UCSC_Control.csv",
    chipsp_file = "ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv",
    out_dir = "Enrichr_dotplots",
    symbol_col = "symbol",
    databases = c(
      "GO_Biological_Process_2023",
      "GO_Cellular_Component_2023",
      "GO_Molecular_Function_2023",
      "KEGG_2021_Human",
      "Reactome_2022",
      "ChEA_2022",
      "MSigDB_Hallmark_2020"
    ),
    top_n_summary = 50L,
    top_n_plot = 10L,
    p_cutoff = 0.05,
    groups_to_plot = c(
      "ChIP",
      "ChIPSP",
      "ChIPSP_only",
      "Shared"
    ),
    uppercase_genes = TRUE,
    seed = 1234L
) {

  # ==========================================================
  # 1. Validate inputs
  # ==========================================================

  if (!is.character(chip_file) || length(chip_file) != 1L ||
      is.na(chip_file) || chip_file == "") {
    stop("`chip_file` must be one valid file path.")
  }

  if (!file.exists(chip_file)) {
    stop("ChIP annotation file not found: ", chip_file)
  }

  if (!is.character(chipsp_file) || length(chipsp_file) != 1L ||
      is.na(chipsp_file) || chipsp_file == "") {
    stop("`chipsp_file` must be one valid file path.")
  }

  if (!file.exists(chipsp_file)) {
    stop("ChIP-SP annotation file not found: ", chipsp_file)
  }

  if (!is.character(out_dir) || length(out_dir) != 1L ||
      is.na(out_dir) || out_dir == "") {
    stop("`out_dir` must be one valid directory path.")
  }

  if (!is.character(symbol_col) || length(symbol_col) != 1L ||
      is.na(symbol_col) || symbol_col == "") {
    stop("`symbol_col` must be one valid column name.")
  }

  if (!is.character(databases) || length(databases) == 0L) {
    stop("`databases` must contain at least one Enrichr database name.")
  }

  if (!is.numeric(top_n_summary) || length(top_n_summary) != 1L ||
      is.na(top_n_summary) || top_n_summary < 1) {
    stop("`top_n_summary` must be one positive integer.")
  }

  if (!is.numeric(top_n_plot) || length(top_n_plot) != 1L ||
      is.na(top_n_plot) || top_n_plot < 1) {
    stop("`top_n_plot` must be one positive integer.")
  }

  if (!is.numeric(p_cutoff) || length(p_cutoff) != 1L ||
      is.na(p_cutoff) || !is.finite(p_cutoff) ||
      p_cutoff <= 0 || p_cutoff > 1) {
    stop("`p_cutoff` must be greater than 0 and less than or equal to 1.")
  }

  if (!is.logical(uppercase_genes) || length(uppercase_genes) != 1L ||
      is.na(uppercase_genes)) {
    stop("`uppercase_genes` must be TRUE or FALSE.")
  }

  # ==========================================================
  # 2. Check required packages
  # ==========================================================

  required_packages <- c(
    "enrichR",
    "ggplot2",
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
      paste(missing_packages, collapse = ", ")
    )
  }

  # ==========================================================
  # 3. Initialize output and Enrichr
  # ==========================================================

  set.seed(seed)

  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  enrichR::setEnrichrSite("Enrichr")

  # ==========================================================
  # 4. Internal helper: read gene list
  # ==========================================================

  get_gene_list <- function(file) {

    tab <- utils::read.csv(
      file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (!symbol_col %in% colnames(tab)) {
      stop(
        "Column `",
        symbol_col,
        "` was not found in: ",
        basename(file)
      )
    }

    genes <- unique(
      as.character(tab[[symbol_col]])
    )

    genes <- trimws(genes)
    genes <- genes[!is.na(genes) & genes != ""]

    if (uppercase_genes) {
      genes <- toupper(genes)
    }

    unique(genes)
  }

  # ==========================================================
  # 5. Internal helper: run Enrichr and save tables
  # ==========================================================

  run_enrichr_save <- function(gene_list, group_name) {

    message("Running Enrichr for: ", group_name)
    message("Input genes: ", length(gene_list))

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
        " has fewer than 5 genes. Skipping."
      )
      return(NULL)
    }

    enr <- tryCatch(
      enrichR::enrichr(
        gene_list,
        databases
      ),
      error = function(e) {
        warning(
          "Enrichr failed for ",
          group_name,
          ": ",
          conditionMessage(e)
        )
        return(NULL)
      }
    )

    if (is.null(enr)) {
      return(NULL)
    }

    for (db in names(enr)) {

      utils::write.csv(
        enr[[db]],
        file = file.path(
          group_dir,
          paste0(
            db,
            "_",
            group_name,
            ".csv"
          )
        ),
        row.names = FALSE,
        quote = FALSE
      )
    }

    enr
  }

  # ==========================================================
  # 6. Internal helper: extract top terms
  # ==========================================================

  extract_top_terms <- function(enr_list, group_name) {

    if (is.null(enr_list)) {
      return(NULL)
    }

    result_list <- lapply(
      names(enr_list),
      function(db) {

        tab <- as.data.frame(
          enr_list[[db]],
          stringsAsFactors = FALSE
        )

        required_columns <- c(
          "Term",
          "P.value",
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
            " for ",
            group_name,
            " because required columns are missing: ",
            paste(missing_columns, collapse = ", ")
          )
          return(NULL)
        }

        tab$Database <- db
        tab$Group <- group_name
        tab$Term_clean <- as.character(tab$Term)

        tab$neg_log10_p <- -log10(
          suppressWarnings(
            as.numeric(tab$P.value)
          )
        )

        if ("Adjusted.P.value" %in% colnames(tab)) {
          tab$neg_log10_adjP <- -log10(
            suppressWarnings(
              as.numeric(tab$Adjusted.P.value)
            )
          )
        } else {
          tab$neg_log10_adjP <- NA_real_
        }

        overlap_text <- as.character(
          tab$Overlap
        )

        tab$overlap_hit <- suppressWarnings(
          as.numeric(
            sub("/.*", "", overlap_text)
          )
        )

        tab$overlap_total <- suppressWarnings(
          as.numeric(
            sub(".*/", "", overlap_text)
          )
        )

        tab$GeneRatio <- tab$overlap_hit /
          tab$overlap_total

        tab <- tab[
          order(
            suppressWarnings(
              as.numeric(tab$P.value)
            ),
            na.last = NA
          ),
          ,
          drop = FALSE
        ]

        if (nrow(tab) > top_n_summary) {
          tab <- tab[
            seq_len(top_n_summary),
            ,
            drop = FALSE
          ]
        }

        tab
      }
    )

    result_list <- Filter(
      Negate(is.null),
      result_list
    )

    if (length(result_list) == 0L) {
      return(NULL)
    }

    out <- do.call(
      rbind,
      result_list
    )

    rownames(out) <- NULL
    out
  }

  # ==========================================================
  # 7. Internal helper: figure theme
  # ==========================================================

  theme_nature_dot <- function(base_size = 8) {

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
  # 8. Internal helper: plot Enrichr dot plot
  # ==========================================================

  plot_enrichr_dotplot <- function(
      summary_all,
      group_name,
      database_name
  ) {

    if (is.null(summary_all) || nrow(summary_all) == 0L) {
      return(NULL)
    }

    p_values <- suppressWarnings(
      as.numeric(summary_all$P.value)
    )

    keep <- (
      summary_all$Group == group_name &
        summary_all$Database == database_name &
        !is.na(p_values) &
        p_values < p_cutoff
    )

    df <- summary_all[
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
      return(NULL)
    }

    df$P.value <- suppressWarnings(
      as.numeric(df$P.value)
    )

    df <- df[
      order(
        df$P.value,
        na.last = NA
      ),
      ,
      drop = FALSE
    ]

    if (nrow(df) > top_n_plot) {
      df <- df[
        seq_len(top_n_plot),
        ,
        drop = FALSE
      ]
    }

    term_clean <- as.character(
      df$Term_clean
    )

    term_clean <- gsub(
      "Homo sapiens",
      "",
      term_clean,
      fixed = TRUE
    )

    term_clean <- gsub(
      "_",
      " ",
      term_clean,
      fixed = TRUE
    )

    term_clean <- gsub(
      "\\s+",
      " ",
      term_clean
    )

    term_clean <- gsub(
      "\\s*R-HSA-[0-9]+",
      "",
      term_clean
    )

    term_clean <- trimws(
      term_clean
    )

    term_clean <- ifelse(
      nchar(term_clean) > 58,
      paste0(
        substr(term_clean, 1, 55),
        "..."
      ),
      term_clean
    )

    df$Term_clean <- term_clean
    df$log10_p <- -log10(df$P.value)

    df$Term_id <- factor(
      seq_len(nrow(df)),
      levels = rev(
        seq_len(nrow(df))
      )
    )

    y_labels <- stats::setNames(
      df$Term_clean,
      as.character(
        seq_len(nrow(df))
      )
    )

    max_ratio <- max(
      df$GeneRatio,
      na.rm = TRUE
    )

    max_log10_p <- max(
      df$log10_p,
      na.rm = TRUE
    )

    if (!is.finite(max_ratio) || max_ratio <= 0) {
      max_ratio <- 1
    }

    if (!is.finite(max_log10_p)) {
      max_log10_p <- 1.01
    }

    color_upper <- max(
      1.01,
      max_log10_p
    )

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = GeneRatio,
        y = Term_id
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
      ggplot2::scale_y_discrete(
        labels = y_labels
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

    pdf_device <- if (isTRUE(capabilities("cairo"))) {
      grDevices::cairo_pdf
    } else {
      "pdf"
    }

    ggplot2::ggsave(
      filename = pdf_file,
      plot = p,
      device = pdf_device,
      width = 6.6,
      height = 3.6,
      units = "in"
    )

    ggplot2::ggsave(
      filename = png_file,
      plot = p,
      width = 6.6,
      height = 3.6,
      units = "in",
      dpi = 600
    )

    p
  }

  # ==========================================================
  # 9. Generate gene sets
  # ==========================================================

  chip_genes <- get_gene_list(
    chip_file
  )

  chipsp_genes <- get_gene_list(
    chipsp_file
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

  gene_sets <- list(
    ChIP = chip_genes,
    ChIPSP = chipsp_genes,
    ChIPSP_only = chipsp_only_genes,
    ChIP_only = chip_only_genes,
    Shared = shared_genes
  )

  # ==========================================================
  # 10. Save gene summary and gene lists
  # ==========================================================

  gene_summary <- data.frame(
    Category = names(gene_sets),
    Gene_count = vapply(
      gene_sets,
      length,
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    gene_summary,
    file = file.path(
      out_dir,
      "Gene_set_summary.csv"
    ),
    row.names = FALSE,
    quote = FALSE
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

  print(gene_summary)

  # ==========================================================
  # 11. Run Enrichr for all groups
  # ==========================================================

  enrichment_results <- list(
    ChIP = run_enrichr_save(
      chip_genes,
      "ChIP"
    ),
    ChIPSP = run_enrichr_save(
      chipsp_genes,
      "ChIPSP"
    ),
    ChIPSP_only = run_enrichr_save(
      chipsp_only_genes,
      "ChIPSP_only"
    ),
    ChIP_only = run_enrichr_save(
      chip_only_genes,
      "ChIP_only"
    ),
    Shared = run_enrichr_save(
      shared_genes,
      "Shared"
    )
  )

  # ==========================================================
  # 12. Merge Enrichr results
  # ==========================================================

  summary_parts <- lapply(
    names(enrichment_results),
    function(group_name) {
      extract_top_terms(
        enrichment_results[[group_name]],
        group_name
      )
    }
  )

  summary_parts <- Filter(
    Negate(is.null),
    summary_parts
  )

  if (length(summary_parts) == 0L) {
    summary_all <- data.frame()
  } else {
    summary_all <- do.call(
      rbind,
      summary_parts
    )
    rownames(summary_all) <- NULL
  }

  utils::write.csv(
    summary_all,
    file = file.path(
      out_dir,
      "Combined_Enrichr_results_all_groups.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )

  # ==========================================================
  # 13. Generate figures
  # ==========================================================

  plot_list <- list()

  if (nrow(summary_all) > 0L) {

    for (group_name in groups_to_plot) {

      for (database_name in databases) {

        plot_key <- paste0(
          group_name,
          "__",
          database_name
        )

        plot_list[[plot_key]] <- plot_enrichr_dotplot(
          summary_all = summary_all,
          group_name = group_name,
          database_name = database_name
        )
      }
    }
  }

  message(
    "Enrichr analysis completed. Output directory: ",
    normalizePath(
      out_dir,
      mustWork = FALSE
    )
  )

  invisible(
    list(
      gene_sets = gene_sets,
      gene_summary = gene_summary,
      enrichment = enrichment_results,
      combined_results = summary_all,
      plots = plot_list
    )
  )
}
