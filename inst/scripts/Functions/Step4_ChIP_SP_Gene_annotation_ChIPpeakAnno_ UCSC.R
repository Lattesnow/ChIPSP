suppressPackageStartupMessages({
  library(ChIPpeakAnno)
  library(GenomicRanges)
  library(dplyr)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(AnnotationDbi)
})

# Decide species and reference genome
BINDING_BP <- 5000

SPECIES <- "human"   # "human" or "mouse"
REF_GENOME <- "hg38" # human: hg19/hg38; mouse: mm9/mm10

stopifnot(SPECIES %in% c("human", "mouse"))

if (SPECIES == "human") {
  stopifnot(REF_GENOME %in% c("hg19", "hg38"))
  
  txdb <- switch(
    REF_GENOME,
    hg19 = TxDb.Hsapiens.UCSC.hg19.knownGene,
    hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene
  )
  
  orgdb <- "org.Hs.eg.db"
  
} else if (SPECIES == "mouse") {
  stopifnot(REF_GENOME %in% c("mm9", "mm10"))
  
  txdb <- switch(
    REF_GENOME,
    mm9  = TxDb.Mmusculus.UCSC.mm9.knownGene,
    mm10 = TxDb.Mmusculus.UCSC.mm10.knownGene
  )
  
  orgdb <- "org.Mm.eg.db"
}

annoData <- toGRanges(txdb, feature = "gene")


# annotate peak table
annotate_peak_table <- function(df, has_score = FALSE, out_csv) {
  req <- c("chr", "start", "end", "pileup")
  if (has_score) req <- c(req, "score")
  if (!all(req %in% colnames(df))) {
    stop("Required columns missing: ", paste(setdiff(req, colnames(df)), collapse = ", "))
  }
  
  gr <- GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start, end = df$end),
    pileup   = df$pileup
  )
  if (has_score) mcols(gr)$score <- df$score
  
  anno <- annotatePeakInBatch(
    gr,
    AnnotationData = annoData,
    output = "both",
    bindingRegion = c(-BINDING_BP, BINDING_BP)
  )
  
  anno <- addGeneIDs(
    anno,
    orgdb,
    IDs2Add = "symbol",
    feature_id_type = "entrez_id"
  )
  
  anno_df <- as.data.frame(anno)
  
  # merge back pileup/score columns
  gr_df <- as.data.frame(gr)[, c("seqnames", "start", "end", "pileup", if (has_score) "score")]
  anno_df <- merge(anno_df, gr_df, by = c("seqnames", "start", "end"), all.x = TRUE)
  

  pileup_col <- intersect(c("pileup", "pileup.x", "pileup.y"), colnames(anno_df))[1]
  if (is.na(pileup_col)) stop("No pileup column found after merge().")
  
  anno_df <- anno_df %>%
    filter(!is.na(symbol) & symbol != "") %>%
    arrange(desc(.data[[pileup_col]])) %>%
    mutate(rank = dense_rank(-.data[[pileup_col]]))
  
  
  write.csv(anno_df, file = out_csv, row.names = FALSE, quote = FALSE)
  message("Wrote: ", out_csv)
  invisible(anno_df)
}

# 1) Conventional ChIP peaks
bed <- read.delim("Combined_ChIP.xls", sep = "\t", stringsAsFactors = FALSE)
annotate_peak_table(
  df = bed,
  has_score = FALSE,
  out_csv = "ChIP_anno_genes_upAdown_UCSC_Control.csv"
)

# 2) ChIP-SP peaks
chipsp <- read.delim("final_ranked_output.xls", sep = "\t", stringsAsFactors = FALSE)
annotate_peak_table(
  df = chipsp,
  has_score = TRUE,
  out_csv = "ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv"
)
