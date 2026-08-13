#Call peak without using MACS2( using threshold method with Bigwig files)
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(data.table)
  library(openxlsx)
})

bw_path <- list.files(
  path = ".",
  pattern = "\\.bw$",
  full.names = TRUE
)

stopifnot(length(bw_path) > 0)
bw_path <- bw_path[1] 

chrom_sizes <- "hg19.chrom.sizes"
out_prefix  <- "hg19"

# Set up this quantile to get solid peaks
top_quantile <- 0.999  
min_peak_bp  <- 100
merge_gap_bp <- 50
score_method <- "max"

# load chrome size file
cs <- fread(chrom_sizes, header = FALSE)
setnames(cs, c("chr", "len"))
keep_chr <- paste0("chr", c(1:22, "X"))
cs <- cs[chr %in% keep_chr]
seqlens <- setNames(cs$len, cs$chr)

gr <- import(bw_path, format="BigWig")
stopifnot(length(gr) > 0)
gr <- gr[is.finite(gr$score)]


lvl <- seqlevels(gr)
if (!any(grepl("^chr", lvl))) {
  lvl2 <- ifelse(grepl("^([0-9]+|X|Y|M)$", lvl), paste0("chr", lvl), lvl)
  seqlevels(gr) <- lvl2
}
seqnames(gr) <- factor(as.character(seqnames(gr)), levels = seqlevels(gr))


gr <- keepSeqlevels(gr, intersect(seqlevels(gr), names(seqlens)), pruning.mode="coarse")
seqlengths(gr) <- seqlens[seqlevels(gr)]


gr <- gr[gr$score > 0]
stopifnot(length(gr) > 0)

message("Intervals after import + rename + >0 filter: ", length(gr))
message("Score range: ", paste(range(gr$score), collapse=" .. "))

# ---- width-weighted quantile threshold
weighted_quantile <- function(x, w, prob) {
  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= prob)[1]]
}

# Optional speed-up for huge GRanges: sample up to 3e6 intervals to estimate threshold
set.seed(1)
n <- length(gr)
idx <- if (n > 3e6) sample.int(n, 3e6) else seq_len(n)

w <- width(gr)[idx]
s <- gr$score[idx]

thr <- weighted_quantile(s, w, top_quantile)
message(sprintf("Threshold (top %.2f%%): %.5f", (1-top_quantile)*100, thr))

# ---- enriched islands + merge
enriched <- gr[gr$score >= thr]
stopifnot(length(enriched) > 0)
message("Enriched intervals: ", length(enriched))

peaks <- reduce(enriched, min.gapwidth = merge_gap_bp + 1L)
peaks <- peaks[width(peaks) >= min_peak_bp]
stopifnot(length(peaks) > 0)
message("Peaks after reduce/filter: ", length(peaks))

# ---- score peaks (max or width-weighted mean pileup)
# IMPORTANT: score against 'enriched' (not full 'gr') for efficiency and consistency
hits <- findOverlaps(peaks, enriched, ignore.strand = TRUE)

dt <- data.table(
  peak_id = queryHits(hits),
  score   = enriched$score[subjectHits(hits)],
  seg_w   = width(pintersect(peaks[queryHits(hits)], enriched[subjectHits(hits)]))
)

if (score_method == "max") {
  peak_score <- dt[, .(pileup = max(score, na.rm = TRUE)), by = peak_id]
} else if (score_method == "mean") {
  peak_score <- dt[, .(pileup = sum(score * seg_w, na.rm = TRUE) / sum(seg_w)), by = peak_id]
} else stop("score_method must be 'max' or 'mean'")

peaks$pileup <- NA_real_
peaks$pileup[peak_score$peak_id] <- peak_score$pileup

# Any peak with no overlap (rare) gets pileup 0
peaks$pileup[is.na(peaks$pileup)] <- 0

# ---- BED score scaling 0..1000
cap <- as.numeric(quantile(peaks$pileup, 0.999, na.rm = TRUE))
if (!is.finite(cap) || cap <= 0) cap <- max(peaks$pileup, na.rm = TRUE)

peaks$bed_score <- as.integer(
  pmax(pmin(round(pmin(peaks$pileup, cap) / cap * 1000), 1000), 0)
)

peaks$name <- sprintf("MYC_pileupPeak_%07d", seq_along(peaks))

# ---- export BED (0-based start)
out_bed <- paste0(out_prefix, ".q99.pileupPeaks.bed")
bed_df <- data.frame(
  chr   = as.character(seqnames(peaks)),
  start = start(peaks) - 1L,
  end   = end(peaks),
  name  = peaks$name,
  score = peaks$bed_score
)
write.table(bed_df, out_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
message("Wrote: ", out_bed)

# export as excel file
out_xls <- paste0(out_prefix, "_q99_pileupPeaks_rawScores_ChIP.xls")

qc_df <- data.frame(
  chr = as.character(seqnames(peaks)),
  start = start(peaks),
  end   = end(peaks),
  length     = width(peaks),
  name         = peaks$name,
  pileup = peaks$pileup,
  bed_score_0_1000 = peaks$bed_score
)

write.table(
  qc_df,
  file = out_xls,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

message("Wrote .xls (tab-delimited): ", out_xls)


