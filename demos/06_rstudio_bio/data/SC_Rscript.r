library(ascend)
library(BiocParallel)
library(gridExtra)

register(MulticoreParam(workers = 4, progressbar=TRUE), default=TRUE)

em_set <- loadCellRanger("data")
em_set

print("Using colInfo...")
colInfo(em_set)

qc_plots <- plotGeneralQC(em_set)
png(filename="data/qcplot1.png")
grid.arrange(qc_plots$libsize_barplot, qc_plots$libsize_histogram, ncol = 1)
dev.off()
png(filename="data/qcplot2.png")
grid.arrange(qc_plots$ngenes_hist, qc_plots$topgenes_violin, ncol = 1)
dev.off()
png(filename="data/qcplot3.png")
grid.arrange(qc_plots$control_hists$Mt, qc_plots$control_hists$Rb, ncol = 2)
dev.off()
png(filename="data/qcplot4.png")
grid.arrange(qc_plots$control_violins$Mt, qc_plots$control_violins$Rb, ncol = 1)
dev.off()

em_set <- normaliseBatches(em_set)

filtered_set <- filterByOutliers(em_set, cell.threshold = 3, control.threshold = 3)
filtered_set <- filterByControl(filtered_set, control = "Mt", pct.threshold = 20)
filtered_set <- filterByControl(filtered_set, control = "Rb", pct.threshold = 50)
filtered_set <- filterLowAbundanceGenes(filtered_set, pct.threshold = 1)

str(progressLog(filtered_set))

filtered_qc_plots <- plotGeneralQC(filtered_set)
png(filename="data/filtered_pc_plots.png")
grid.arrange(filtered_qc_plots$libsize_histogram, filtered_qc_plots$ngenes_hist, filtered_qc_plots$control_hists$Mt, filtered_qc_plots$control_hists$Rb, filtered_qc_plots$control_violins$Mt, filtered_qc_plots$control_violins$Rb, ncol = 2)
dev.off()

norm_set <- normaliseByRLE(filtered_set)

counts(norm_set)[1:5,1:5]
normcounts(norm_set)[1:5,1:5]

norm_qc <- plotNormQC(norm_set)
png(filename="data/norm_qcplot1.png")
grid.arrange(norm_qc$libsize_histograms$count, norm_qc$libsize_histograms$normcount, ncol = 1)
dev.off()
png(filename="data/norm_qcplot2.png")
grid.arrange(norm_qc$sampled_genes$GAPDH$counts, norm_qc$sampled_genes$GAPDH$normcounts, ncol = 1)
dev.off()
png(filename="data/norm_qcplot3.png")
grid.arrange(norm_qc$sampled_genes$MALAT1$counts, norm_qc$sampled_genes$MALAT1$normcounts, ncol = 1)
dev.off()
png(filename="data/norm_qcplot4.png")
grid.arrange(norm_qc$sampled_cell_gene_expression$counts, norm_qc$sampled_cell_gene_expression$normcounts, ncol = 1)
dev.off()
png(filename="data/filtered_qc_boxplot.png")
filtered_qc_plots$topgenes_boxplot
dev.off()

norm_set <- excludeControl(norm_set, control = c("Mt", "Rb"))
print(plotTopGenesBoxplot(norm_set, n = 20))

pca_set <- runPCA(norm_set, ngenes = 1500, scaling = TRUE)

reducedDim(pca_set, "PCA")[1:5,1:5]

png(filename="data/pca_var_plot.png")
plotPCAVariance(pca_set, n = 50)
dev.off()
png(filename="data/pca_plot.png")
plotPCA(pca_set, PCX = 1, PCY = 2, group = "batch")
dev.off()
clustered_set<-runCORE(pca_set, conservative = FALSE, nres=40, remove.outlier = FALSE)
cluster_analysis <- clusterAnalysis(clustered_set)
cluster_analysis$keyStats

png(filename="data/stability_dendro_plot.png")
plotStabilityDendro(clustered_set)
dev.off()
png(filename="data/stability_plot.png")
plotStability(clustered_set)
dev.off()
png(filename="data/dendogram.png")
plotDendrogram(clustered_set)
dev.off()

clustered_set <- runTSNE(clustered_set, PCA = TRUE, dimensions = 2, seed = 1, perplexity = 30, theta = 0.5)
tsne_plot <- plotTSNE(clustered_set, Dim1 = 1, Dim2 = 2, group = "cluster")
pca_plot <- plotPCA(clustered_set, PCX=1, PCY=2, group = "cluster")
mds_plot <- plotMDS(clustered_set, Dim1 = 1, Dim2 = 2, group = "cluster")
png(filename="data/tsne_plot.png")
tsne_plot
dev.off()
png(filename="data/pca_plot.png")
pca_plot
dev.off()
png(filename="data/mds_plot.png")
mds_plot
dev.off()

col_info <- colInfo(clustered_set)
print("Cells per batch")
table(col_info$batch)
print("Cells per cluster")
table(col_info$cluster)
cluster1_vs_others <- runDiffExpression(clustered_set, group = "cluster", condition.a = 1, condition.b = c(2, 3), subsampling = FALSE, ngenes = 1500)
cluster2_vs_others <- runDiffExpression(clustered_set, group = "cluster", condition.a = 2, condition.b = c(1, 3), subsampling = FALSE, ngenes = 1500)
cluster3_vs_others <- runDiffExpression(clustered_set, group = "cluster", condition.a = 3, condition.b = c(1, 2), subsampling = FALSE, ngenes = 1500)
cluster1_volcano <- plotVolcano(cluster1_vs_others, l2fc = 2, threshold = 5e-3, labels = TRUE, label.size = 3, check.overlap = TRUE)
cluster2_volcano <- plotVolcano(cluster2_vs_others, l2fc = 2, threshold = 5e-3, labels = TRUE, label.size = 3, check.overlap = TRUE)
cluster3_volcano <- plotVolcano(cluster3_vs_others, l2fc = 2, threshold = 5e-3, labels = TRUE, label.size = 3, check.overlap = TRUE)
library(ggplot2)
cluster1_volcano <- cluster1_volcano + ggtitle("Cluster 1 vs Cluster 2 and 3")
cluster2_volcano <- cluster2_volcano + ggtitle("Cluster 2 vs Cluster 1 and 3")
cluster3_volcano <- cluster3_volcano + ggtitle("Cluster 3 vs Cluster 1 and 2")
png(filename="data/cluster1_volcano.png")
cluster1_volcano
dev.off()
png(filename="data/cluster2_volcano.png")
cluster2_volcano
dev.off()
png(filename="data/cluster3_volcano.png")
cluster3_volcano
dev.off()

sce_object <- EMSet2SCE(clustered_set)
count_matrix <- counts(clustered_set)
normalised_matrix <- normcounts(clustered_set)
log_matrix <- logcounts(clustered_set)
scran_normalised <- scranNormalise(filtered_set, quickCluster = FALSE, min_mean = 1e-05)
cluster2_vs_others <- runDESeq(clustered_set, group = "cluster", condition.a = 2, condition.b = c(1, 3), ngenes = 1500)
png(filename="data/volcano_deseq.png")
plotVolcano(cluster1_vs_others, labels = TRUE)
dev.off()
cluster1_vs_others <- runDESeq2(clustered_set, group = "cluster", condition.a = 1, condition.b = c(2, 3), ngenes = 1500)
