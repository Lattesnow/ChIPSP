# ============================================================
# Enrichr analysis dotplots
# ChIP, ChIP-SP, and ChIP-SP-only genes
# ============================================================

library(enrichR)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(Cairo)

set.seed(1234)

# inputs from previous function
chip_file   <- "ChIP_anno_genes_upAdown_UCSC_Control.csv"
chipsp_file <- "ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv"

out_dir <- "Enrichr_dotplots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# enrichdatabase call
setEnrichrSite("Enrichr")

dbs <- c(
  "GO_Biological_Process_2023",
  "GO_Cellular_Component_2023",
  "GO_Molecular_Function_2023",
  "KEGG_2021_Human",
  "Reactome_2022",
  "ChEA_2022",
  "MSigDB_Hallmark_2020"
)

# ============================================================
# Functions
# ============================================================

get_gene_list <- function(file, symbol_col = "symbol") {
  
  tab <- read.csv(file, stringsAsFactors = FALSE)
  
  genes <- unique(tab[[symbol_col]])
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- toupper(genes)
  
  return(genes)
}


run_enrichr_save <- function(gene_list, group_name, dbs, out_dir) {
  
  message("Running Enrichr for: ", group_name)
  message("Input genes: ", length(gene_list))
  
  group_dir <- file.path(out_dir, group_name)
  
  dir.create(
    group_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  if (length(gene_list) < 5) {
    
    warning(
      group_name,
      " has fewer than 5 genes. Skipping."
    )
    
    return(NULL)
  }
  
  enr <- enrichr(gene_list, dbs)
  
  for (db in names(enr)) {
    
    write_csv(
      enr[[db]],
      file.path(
        group_dir,
        paste0(
          db,
          "_",
          group_name,
          ".csv"
        )
      )
    )
  }
  
  return(enr)
}


extract_top_terms <- function(enr_list, group_name, top_n = 50) {
  
  if (is.null(enr_list)) return(NULL)
  
  bind_rows(
    lapply(
      names(enr_list),
      function(db) {
        
        enr_list[[db]] %>%
          as_tibble() %>%
          mutate(
            Database = db,
            Group = group_name,
            Term_clean = Term,
            
            neg_log10_p = -log10(P.value),
            neg_log10_adjP = -log10(Adjusted.P.value),
            
            # Number of input genes overlapping the pathway
            overlap_hit = as.numeric(
              sub("/.*", "", Overlap)
            ),
            
            # Total number of genes in the pathway
            overlap_total = as.numeric(
              sub(".*/", "", Overlap)
            ),
            
            # Gene ratio:
            # overlapping genes / total pathway genes
            GeneRatio = overlap_hit / overlap_total
          ) %>%
          
          # Rank by P value
          arrange(P.value) %>%
          
          slice_head(n = top_n)
      }
    )
  )
}


# ============================================================
# Arial font
# ============================================================

theme_nature_dot <- function(base_size = 8) {
  
  theme_classic(
    base_size = base_size,
    base_family = "Arial"
  ) +
    theme(
      
      text = element_text(
        family = "Arial",
        color = "black"
      ),
      
      plot.title = element_text(
        family = "Arial",
        size = 10,
        face = "bold",
        hjust = 0
      ),
      
      axis.title = element_text(
        family = "Arial",
        size = 9
      ),
      
      axis.text = element_text(
        family = "Arial",
        size = 8,
        color = "black"
      ),
      
      axis.text.x = element_text(
        family = "Arial",
        size = 8,
        angle = 45,
        hjust = 1
      ),
      
      axis.text.y = element_text(
        family = "Arial",
        size = 10,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.5,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.5,
        color = "black"
      ),
      
      panel.grid.major = element_line(
        linewidth = 0.25,
        color = "grey88"
      ),
      
      panel.grid.minor = element_blank(),
      
      legend.title = element_text(
        family = "Arial",
        size = 8
      ),
      
      legend.text = element_text(
        family = "Arial",
        size = 7
      ),
      
      plot.margin = margin(
        5,
        8,
        5,
        5
      )
    )
}


# ============================================================
# Plot Enrichr dotplot
#
# Pathway ranking = P value
# Dot color       = -log10(P value)
# Dot size        = gene count
# X-axis          = overlap_hit / overlap_total
# ============================================================

plot_enrichr_dotplot <- function(
    summary_all,
    group_name,
    database_name,
    top_n = 10,
    out_dir,
    p_cutoff = 0.05
) {
  
  df <- summary_all %>%
    
    filter(
      Group == group_name,
      Database == database_name
    ) %>%
    
    filter(
      P.value < p_cutoff
    ) %>%
    
    # --------------------------------------------------------
  # Rank pathways by P value
  # smallest P value first
  # --------------------------------------------------------
  arrange(P.value) %>%
    
    slice_head(n = top_n) %>%
    
    mutate(
      
      Term_clean = str_replace(
        Term_clean,
        "Homo sapiens",
        ""
      ),
      
      Term_clean = str_replace_all(
        Term_clean,
        "_",
        " "
      ),
      
      Term_clean = str_replace_all(
        Term_clean,
        "\\s+",
        " "
      ),
      
      Term_clean = str_replace(
        Term_clean,
        "\\s*R-HSA-[0-9]+",
        ""
      ),
      
      Term_clean = str_trim(Term_clean),
      
      Term_clean = str_trunc(
        Term_clean,
        58
      ),
      
      # ------------------------------------------------------
      # Put smallest P-value pathway at top
      # ------------------------------------------------------
      Term_clean = factor(
        Term_clean,
        levels = rev(Term_clean)
      ),
      
      # ------------------------------------------------------
      # Dot color
      # ------------------------------------------------------
      log10_p = -log10(P.value)
    )
  
  
  if (nrow(df) == 0) {
    
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
  
  
  max_ratio <- max(
    df$GeneRatio,
    na.rm = TRUE
  )
  
  max_log10_p <- max(
    df$log10_p,
    na.rm = TRUE
  )
  
  
  # ==========================================================
  # Dotplot
  # ==========================================================
  
  p <- ggplot(
    df,
    aes(
      x = GeneRatio,
      y = Term_clean
    )
  ) +
    
    # --------------------------------------------------------
  # Dot size = gene count
  # Dot color = -log10(P)
  # --------------------------------------------------------
  geom_point(
    aes(
      size = overlap_hit,
      color = log10_p
    ),
    alpha = 0.95
  ) +
    
    
    # --------------------------------------------------------
  # P-value color
  # --------------------------------------------------------
  scale_color_gradientn(
    colors = c(
      "#2F8DD8",
      "#9B8FBF",
      "#E46A6A"
    ),
    
    limits = c(
      1,
      max_log10_p
    ),
    
    oob = scales::squish,
    
    name = expression(
      -log[10](P)
    )
  ) +
    
    
    # --------------------------------------------------------
  # Dot size = number of overlapping genes
  # NOT odds ratio
  # --------------------------------------------------------
  scale_size_continuous(
    range = c(
      1.8,
      6.5
    ),
    
    breaks = scales::pretty_breaks(
      n = 4
    ),
    
    name = "Gene count"
  ) +
    
    
    # --------------------------------------------------------
  # X-axis
  #
  # GeneRatio =
  # overlapping genes / total pathway genes
  # --------------------------------------------------------
  scale_x_continuous(
    limits = c(
      0,
      max_ratio * 1.12
    ),
    
    expand = expansion(
      mult = c(
        0.03,
        0.08
      )
    )
  ) +
    
    
    labs(
      title = paste0(
        group_name,
        " - ",
        database_name
      ),
      
      x = "Gene ratio",
      
      y = NULL
    ) +
    
    
    theme_nature_dot() +
    
    
    theme(
      plot.margin = margin(
        5,
        35,
        5,
        20
      ),
      
      legend.position = "right"
    )
  
  
  # ==========================================================
  # Output
  # ==========================================================
  
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
  
  
  ggsave(
    file.path(
      out_dir,
      paste0(
        "dotplot_",
        safe_group,
        "_",
        safe_db,
        ".pdf"
      )
    ),
    
    p,
    
    device = cairo_pdf,
    
    width = 6.6,
    height = 3.6,
    units = "in"
  )
  
  
  ggsave(
    file.path(
      out_dir,
      paste0(
        "dotplot_",
        safe_group,
        "_",
        safe_db,
        ".png"
      )
    ),
    
    p,
    
    width = 6.6,
    height = 3.6,
    units = "in",
    dpi = 600
  )
  
  
  return(p)
}


# ============================================================
# Input gene lists
# ============================================================

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


# ============================================================
# Gene summary
# ============================================================

gene_summary <- tibble(
  
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
  )
)


write_csv(
  gene_summary,
  file.path(
    out_dir,
    "Gene_set_summary.csv"
  )
)


write_lines(
  chip_genes,
  file.path(
    out_dir,
    "ChIP_gene_list.txt"
  )
)


write_lines(
  chipsp_genes,
  file.path(
    out_dir,
    "ChIPSP_gene_list.txt"
  )
)


write_lines(
  chipsp_only_genes,
  file.path(
    out_dir,
    "ChIPSP_only_gene_list.txt"
  )
)


write_lines(
  chip_only_genes,
  file.path(
    out_dir,
    "ChIP_only_gene_list.txt"
  )
)


write_lines(
  shared_genes,
  file.path(
    out_dir,
    "Shared_gene_list.txt"
  )
)


print(gene_summary)


# ============================================================
# EnrichR
# ============================================================

enr_chip <- run_enrichr_save(
  chip_genes,
  "ChIP",
  dbs,
  out_dir
)


enr_chipsp <- run_enrichr_save(
  chipsp_genes,
  "ChIPSP",
  dbs,
  out_dir
)


enr_chipsp_only <- run_enrichr_save(
  chipsp_only_genes,
  "ChIPSP_only",
  dbs,
  out_dir
)


enr_chip_only <- run_enrichr_save(
  chip_only_genes,
  "ChIP_only",
  dbs,
  out_dir
)


enr_shared <- run_enrichr_save(
  shared_genes,
  "Shared",
  dbs,
  out_dir
)


# ============================================================
# Merge
# ============================================================

summary_all <- bind_rows(
  
  extract_top_terms(
    enr_chip,
    "ChIP"
  ),
  
  extract_top_terms(
    enr_chipsp,
    "ChIPSP"
  ),
  
  extract_top_terms(
    enr_chipsp_only,
    "ChIPSP_only"
  ),
  
  extract_top_terms(
    enr_chip_only,
    "ChIP_only"
  ),
  
  extract_top_terms(
    enr_shared,
    "Shared"
  )
)


write_csv(
  summary_all,
  file.path(
    out_dir,
    "Combined_Enrichr_results_all_groups.csv"
  )
)


# ============================================================
# Generate figures
# ============================================================

groups_to_plot <- c(
  "ChIP",
  "ChIPSP",
  "ChIPSP_only",
  "Shared"
)

databases_to_plot <- dbs


for (g in groups_to_plot) {
  
  for (db in databases_to_plot) {
    
    plot_enrichr_dotplot(
      
      summary_all = summary_all,
      
      group_name = g,
      
      database_name = db,
      
      top_n = 10,
      
      out_dir = out_dir,
      
      p_cutoff = 0.05
    )
  }
}