library(ClusterGVis)
library(readxl)
library(biomaRt)
library(dplyr)
library(data.table)
# load data
rna_file <- list.files(
  path = ".",
  pattern = "RNA.*\\.(txt|xls|tsv|csv)$",
  ignore.case = TRUE,
  full.names = TRUE
)

rna_file

stopifnot(length(rna_file) == 1)
rna_file <- rna_file[1]

expr <- if (grepl("\\.csv$", rna_file, ignore.case = TRUE)) {
  read.csv(rna_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  read.delim(rna_file, stringsAsFactors = FALSE, check.names = FALSE)
}

expr <- expr[rowSums(expr[, -1, drop = FALSE]) != 0, ]
colnames(expr)[1]<- "gene_name"

exps <- expr %>% distinct(gene_name, .keep_all = TRUE)
rownames(exps) <- exps$gene_name
exps <- exps[,-1]

total_reads_per_sample <- colSums(exps)
exps <- sweep(exps, 2, total_reads_per_sample, "/") * 1e6  

write.csv(exps, file = "expression_normalized_RPM.csv")

tab_ChIP   <- fread("ChIP_anno_genes_upAdown_UCSC_Control.csv")
tab_CHIPSP <- fread("ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv")

gene_list_ChIP     <- unique(na.omit(tab_ChIP$symbol))
gene_list_ChIP_SP  <- unique(na.omit(tab_CHIPSP$symbol))

ChIP_SP_Shared <- sort(intersect(gene_list_ChIP_SP, gene_list_ChIP))
ChIP_all       <- sort(gene_list_ChIP)
ChIP_SP_all    <- sort(gene_list_ChIP_SP)
ChIP_SP_Unique <- sort(setdiff(gene_list_ChIP_SP, gene_list_ChIP))

# Pad to equal length for "four columns"
pad_to <- function(x, n) c(x, rep(NA_character_, n - length(x)))
max_n <- max(length(ChIP_SP_Shared), length(ChIP_all), length(ChIP_SP_all), length(ChIP_SP_Unique))

out_df <- data.frame(
  ChIP_SP_Shared = pad_to(ChIP_SP_Shared, max_n),
  ChIP           = pad_to(ChIP_all, max_n),
  ChIP_SP        = pad_to(ChIP_SP_all, max_n),
  ChIP_SP_Unique = pad_to(ChIP_SP_Unique, max_n),
  stringsAsFactors = FALSE
)

# Optional: write out
write.table(out_df, file = "ChIP_vs_ChIPSP_geneSetOverlap_4columns.xls",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Quick counts
cat("Counts:\n",
    "ChIP =", length(ChIP_all), "\n",
    "ChIP_SP =", length(ChIP_SP_all), "\n",
    "Shared =", length(ChIP_SP_Shared), "\n",
    "ChIP_SP_Unique =", length(ChIP_SP_Unique), "\n")

filtered_results <- list()

ref <- out_df

exp_df <- as.data.frame(exps)
exp_df$gene_id <- rownames(exp_df)

gene_sets <- c("ChIP_SP_Shared","ChIP", "ChIP_SP", "ChIP_SP_Unique" )

# define top genes
get_top30 <- function(tab) {
  tab <- as.data.frame(tab)
  tab$symbol <- as.character(tab$symbol)
  tab <- tab[!is.na(tab$symbol) & tab$symbol != "", ]
  
  if ("rank" %in% names(tab)) {
    tab <- tab[order(tab$rank, decreasing = FALSE), ]
  } else if ("pileup" %in% names(tab)) {
    tab <- tab[order(tab$pileup, decreasing = TRUE), ]
  } else if ("pileup.x" %in% names(tab)) {
    tab <- tab[order(tab$pileup.x, decreasing = TRUE), ]
  }
  
  unique(tab$symbol)[1:min(30, length(unique(tab$symbol)))]
}

top30_CHIP   <- get_top30(tab_ChIP)
top30_CHIPSP <- get_top30(tab_CHIPSP)

top30_SHARED <- intersect(top30_CHIP, top30_CHIPSP)  # per your rule

out_dir <- "ClusterGVis_heatmaps"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (gene_set in gene_sets) {

# Convert rownames to a column to use dplyr functions


  exp_filtered_df <- exp_df %>%
    filter(gene_id %in% ref[[gene_set]])

rownames(exp_filtered_df) <- exp_filtered_df$gene_id
exp_filtered_df$gene_id <- NULL 

exp_filtered <- as.matrix(exp_filtered_df)
exps <- as.data.frame(exp_filtered)

# exps <- exps[rowSums(exps) > 200, ]

# check optimal cluster numbers
getClusters(exps)

# using mfuzz for clustering
cm <- clusterData(exps,
                  cluster.method = "mfuzz",
                  cluster.num = 6)
# using TCseq for clustering
ct <- clusterData(exps,
                  cluster.method = "TCseq",
                  cluster.num = 6)

# using kemans for clustering
ck <- clusterData(exps,
                  cluster.method = "kmeans",
                  cluster.num = 6)


# supply other aruguments passed by Heatmap function
visCluster(object = ck,
           plot.type = "heatmap",
           column_names_rot = 45)



# Choose mark list per your rules, then intersect with current matrix rownames
if (gene_set %in% c("ChIP_SP", "ChIP_SP_Unique")) {
  mark_source <- top30_CHIPSP
} else if (gene_set == "ChIP") {
  mark_source <- top30_CHIP
} else if (gene_set == "ChIP_SP_Shared") {
  mark_source <- top30_SHARED
} else {
  mark_source <- character(0)
}

markGenes <- intersect(rownames(exps), mark_source)
markGenes <- markGenes[1:min(30, length(markGenes))]  # safe cap

pdf(file.path(out_dir, paste0("addgene_CM_", gene_set, ".pdf")),
    height = 8, width = 6, onefile = FALSE)
visCluster(object = cm,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           )
dev.off()


pdf(file.path(out_dir, paste0("addgene_CK_", gene_set, ".pdf")),
    height = 8, width = 6, onefile = FALSE)
visCluster(object = ck,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
          )
dev.off()
}