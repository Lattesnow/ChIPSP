parent_dir <- getwd()
topdirs <- list.dirs(parent_dir, recursive = FALSE, full.names = TRUE)

count_unique_genes <- function(target_file) {
  
  res <- lapply(topdirs, function(td) {
    
    files <- list.files(
      td,
      pattern = paste0("^", target_file, "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(files) == 0) {
      return(c(
        n_files_found = 0L,
        total_rows = NA_integer_,
        n_unique_genes = NA_integer_
      ))
    }
    
    all_symbols <- unlist(lapply(files, function(f) {
      df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df) || !"symbol" %in% colnames(df)) return(character(0))
      syms <- df$symbol
      syms[!is.na(syms) & syms != ""]
    }))
    
    total_rows <- sum(sapply(files, function(f) {
      df <- tryCatch(read.csv(f), error = function(e) NULL)
      if (is.null(df)) return(0L)
      nrow(df)
    }))
    
    c(
      n_files_found = length(files),
      total_rows = total_rows,
      n_unique_genes = length(unique(all_symbols))
    )
  })
  
  res <- as.data.frame(do.call(rbind, res))
  rownames(res) <- basename(topdirs)
  res
}

# Run two output files
chipsp_res  <- count_unique_genes("ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv")
control_res <- count_unique_genes("ChIP_anno_genes_upAdown_UCSC_Control.csv")

colnames(chipsp_res)  <- paste0("CHIPSP_", colnames(chipsp_res))
colnames(control_res) <- paste0("Control_", colnames(control_res))

final_res <- cbind(chipsp_res, control_res)

# Include the HiC and ChIP files rows here

count_rows_for_file <- function(target_file) {
  
  res <- lapply(topdirs, function(td) {
    
    files <- list.files(
      td,
      pattern = paste0("^", target_file, "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(files) == 0) {
      return(NA_integer_)
    }
    
    total_rows <- sum(sapply(files, function(f) {
      df <- tryCatch(read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df)) return(0L)
      nrow(df)
    }))
    
    total_rows
  })
  
  res <- unlist(res)
  names(res) <- basename(topdirs)
  res
}

# Count rows
chip_rows <- count_rows_for_file("Combined_ChIP.xls")
hic_rows  <- count_rows_for_file("Combined_HiC.xls")

# Add columns
final_res$Combined_ChIP_rows <- chip_rows[rownames(final_res)]
final_res$Combined_HiC_rows  <- hic_rows[rownames(final_res)]


out_file <- "CHIPSP_vs_Control_unique_gene_counts.csv"
write.csv(final_res, out_file, row.names = TRUE)

print(final_res)
cat("Wrote:", out_file, "\n")
