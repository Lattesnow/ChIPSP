## Load dataset and make metadata
counts = read.delim("readcount.xls", header=1,row.names = 1)
#counts2 <- counts[,-1]
#rownames(counts2) <- counts[,1]

library(DESeq2)


meta = sapply(colnames(counts), function (x) substring(x,1,1))
meta = data.frame(meta)
colnames(meta) = c("Condition")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design= ~ Condition)

filt <- rowSums(counts(dds) < 10) > dim(meta)[1]*0.9
dds <- dds[!filt,]

## Perform DESeq2()
dds = DESeq(dds)
res = results(dds)
write.table(res,file = "Deseq2_output.txt",sep = "\t")


library("org.Hs.eg.db")

df1<-read.delim("out_count.xls",sep = "\t",row.names = 1)
df1<- res

data.df<- as.data.frame(df1)
data.df$symbol<- mapIds(org.Hs.eg.db,keys = rownames(data.df),keytype = "ENSEMBL",column = "SYMBOL")
write.table(data.df,file = "output_SYMBOL.txt",sep = "\t")