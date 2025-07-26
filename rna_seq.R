count_files <-c("kidney_a.s.count.txt","kidney_b.s.count.txt","kidney_c.s.count.txt","skin_5e.s.count.txt","skin_5f.s.count.txt","skin_6a.s.count.txt")
count_list <- lapply(count_files, function(file) {read.table(file, header = TRUE, row.names = 1, sep = "\t")})
counts <- do.call(cbind, count_list)
colnames(counts) <- c("kidney_a","kidney_b","kidney_c","skin_5e", "skin_5f","skin_6a")
counts <- counts[, c("kidney_a", "kidney_b", "kidney_c", "skin_5e", "skin_5f", "skin_6a")]
colnames(counts)
head(counts)
group <- factor(c("kidney", "kidney", "kidney","skin", "skin", "skin"))
y <- DGEList(counts = counts, group = group)
y
y <- calcNormFactors(y)
plotMDS(y)
y <- estimateDisp(y)
plotBCV(y)
et <- exactTest(y, pair = c("kidney", "skin"))
de <- decideTests(et)
summary(de)
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=0, col="blue")
diffExpGenes <- topTags(et, n=1000, p.value = 0.05)
head(diffExpGenes$table)
topTags(et)
sample_data <- data.frame(
  sample_name = c("kidney_a", "kidney_b", "kidney_c", "skin_5e", "skin_5f", "skin_6a"),
  Group = factor(c("kidney", "kidney", "kidney", "skin", "skin", "skin"))
)
rownames(sample_data) <- sample_data$sample_name
data_deseq <- DESeqDataSetFromMatrix(countData = counts, colData = sample_data, design = ~ Group)
data_deseq
nrow(data_deseq)
data_deseq <- data_deseq[ rowSums(counts(data_deseq)) > 1, ]
data_deseq
rld <- rlog(data_deseq, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( rld$cell_type, rld$dev_stage, rld$replicate, sep="-" )
sampleLabels <- colnames(rld)
rownames(sampleDistMatrix) <- sampleLabels
colnames(sampleDistMatrix) <- sampleLabels
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col=colors)
plotPCA(rld,intgroup = "Group")
dds <- DESeq(data_deseq)
res <- results(dds)
res <- res[order(res$padj), ]
res <- res[!is.na(res$padj), ]
top50_genes <- rownames(res)[1:50]
rld <- rlog(dds, blind=FALSE)
head(top50_genes)
pheatmap(assay(rld)[top50_genes, ], show_rownames = TRUE, fontsize_row = 7, fontsize_col = 10, scale = "row")
heatmap_matrix <- assay(rld)[topVarGenes,]
annotation_col <- as.data.frame(colData(rld)[,"Group",drop =FALSE])
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)
heatmap_matrix <- assay(rld)[topVarGenes, ]
annotation_col <- as.data.frame(colData(rld)[, "Group", drop = FALSE])
colnames(annotation_col) <- "Tissue_Type"
pheatmap(
  heatmap_matrix,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 12,
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(255),
  main = "Top 50 DE Genes: Kidney vs Skin"
)