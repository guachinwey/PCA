library(Seurat)
library(harmony)
library(Matrix)
library(Rcpp)
library(dplyr)
library(patchwork)
library(BiocManager)
library(stringr)
library(glmGamPoi) # Can contribute to a more efficient workflow.
library(ggplot2)
library(writexl)
library(readxl)
library(pheatmap)
library(ComplexHeatmap)
library(openxlsx)
library(HGNChelper)

merged.datasets <- readRDS("/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat_3.rds")
merged.datasets.markers <- read_excel("/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/FindAllMarkers.xlsx")

top2 <- merged.datasets.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write_xlsx(top2, "/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/TOP2.xlsx")
#top2 <- read_excel("/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/TOP2.genes.xlsx")

top10 <- merged.datasets.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(merged.datasets, features = top10$gene) + NoLegend()

write_xlsx(top10, "/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/TOP10.xlsx")
#top10 <- read_excel("/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/Top10.xlsx")
pdf("/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/heatmap.pdf")

DoHeatmap(merged.datasets, features = top10$gene) + NoLegend()
dev.off()

# As heatmap is too dense, divide by clusters to facilitate the visualization: 
pdf("/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/heatmap.by.clusters.pdf", width = 12, height = 12)

identidades <- levels(Idents(merged.datasets))

cluster_1 <- identidades[0:8]
cluster_2 <- identidades[9:18]
cluster_3 <- identidades[19:27]
cluster_4 <- identidades[28:37]
cluster_5 <- identidades[38:44]

genes_group1 <- top10$gene[top10$cluster %in% cluster_1]
genes_group2 <- top10$gene[top10$cluster %in% cluster_2]
genes_group3 <- top10$gene[top10$cluster %in% cluster_3]
genes_group4 <- top10$gene[top10$cluster %in% cluster_4]
genes_group5 <- top10$gene[top10$cluster %in% cluster_5]

heatmap1 <- DoHeatmap(merged.datasets, group.by = "ident", features = genes_group1, cells = WhichCells(merged.datasets, idents = cluster_1)) + NoLegend() 
heatmap2 <- DoHeatmap(merged.datasets, group.by = "ident", features = genes_group2, cells = WhichCells(merged.datasets, idents = cluster_2)) + NoLegend() 
heatmap3 <- DoHeatmap(merged.datasets, group.by = "ident", features = genes_group3, cells = WhichCells(merged.datasets, idents = cluster_3)) + NoLegend()
heatmap4 <- DoHeatmap(merged.datasets, group.by = "ident", features = genes_group4, cells = WhichCells(merged.datasets, idents = cluster_4)) + NoLegend() 
heatmap5 <- DoHeatmap(merged.datasets, group.by = "ident", features = genes_group5, cells = WhichCells(merged.datasets, idents = cluster_5)) + NoLegend()

print(heatmap1)
print(heatmap2)
print(heatmap3)
print(heatmap4)
print(heatmap5)

dev.off()

#List of genes by population.
db_ <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/characteristic genes.xlsx";
dab <- openxlsx::read.xlsx("/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/characteristic genes.xlsx")
tissue <- "Prostate"
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

gs_list <- gene_sets_prepare(db_, tissue)

es.max <- sctype_score(scRNAseqData = merged.datasets[["SCT"]]@scale.data, scaled = TRUE,
                       gs = gs_list$gs_positive, gs2 = NULL)

cL_results <- do.call("rbind", lapply(unique(merged.datasets@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ , rownames(merged.datasets@meta.data[merged.datasets@meta.data$seurat_clusters==cl, ])]), 
                   decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl,
                  ncells = sum(merged.datasets@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_results %>% 
  group_by(cluster) %>%
  top_n(n = 1, wt = scores)

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

tablaclusters <- sctype_scores[, 1:3]
View(tablaclusters)

write_xlsx(tablaclusters, '/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/AnnotatedClusters.xlsx')

merged.datasets@meta.data$customclassif = ""
for (j in unique(sctype_scores$cluster)) {
  cl_type = sctype_scores[sctype_scores$cluster == j,];
  merged.datasets@meta.data$customclassif[merged.datasets@meta.data$seurat_clusters == j] = 
    as.character(cl_type$type[1])
}

colnames(merged.datasets@meta.data)[colnames(merged.datasets@meta.data) == "customclassif"] <- "clusters"

#Visualization
pdf("/home/csic/epi/sbb/Silvia/Analysis/NEWS/GRAPHICS/umap_population.pdf", width = 12, height = 12)

p1 <- DimPlot(merged.datasets, reduction = "umap", label = TRUE, repel = TRUE, group.by = "clusters")
Idents(object = merged.datasets) <- "sample"
p1
p2 <- DimPlot(merged.datasets, reduction = "umap", repel = TRUE)
p2
p1+p2

dev.off()


pdf("/home/csic/epi/sbb/Silvia/Analysis/NEWS/GRAPHICS/umap_population_by_sample.pdf", width = 12, height = 12)

samples <- unique(merged.datasets$sample)
for (sample in samples) {
  p <- DimPlot(merged.datasets, reduction = "umap", label = TRUE, group.by = "clusters", cells = WhichCells(merged.datasets, expression = sample == !!sample)) +
    ggtitle(paste("Sample:", sample))+
    theme_minimal()
  print(p)
}

dev.off()