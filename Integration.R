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

# Load the Seurat objects
input_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/SAMPLES/"
file_list <- list.files(path = input_dir, pattern = "*.rds", full.names = TRUE)
object_names <- c("CRPC_1", "CRPC_2", "CRPC_3", "CRPC_4", "CRPC_5", "CRPC_6", 
                  "HEALTH1", "HEALTH2", "HEALTH3", "HEALTH4", "HEALTH5", "HEALTH6", "PCA")

seurat_objects <- lapply(file_list, readRDS)
names(seurat_objects) <- object_names
for (i in 1:13) {
  levels(seurat_objects[[i]]@active.ident) <- object_names[i]
}
list2env(seurat_objects, .GlobalEnv)

length(seurat_objects)
lapply(seurat_objects, function(x) levels(x@active.ident))
lapply(seurat_objects, function(x) dim(x))

# Merge all datasets
merged.datasets <- merge(seurat_objects[[1]], y = seurat_objects[-1], project = "PCa_Human")
head(colnames(merged.datasets))
head(merged.datasets)
head(rownames(merged.datasets))
table(merged.datasets@active.ident)
ncol(merged.datasets)

# Pre-processing
merged.datasets[["percent.mt"]] <- PercentageFeatureSet(merged.datasets, pattern = "^MT-")

# Violin plot
VlnPlot(merged.datasets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filtering
minCov <- 1000
countLOW <- ifelse(min(merged.datasets$nCount_RNA) >= minCov, min(merged.datasets$nCount_RNA), quantile(merged.datasets$nCount_RNA, prob = c(0.01)))
countHIGH <- quantile(merged.datasets$nCount_RNA, prob = 0.99)
featureLOW <- quantile(merged.datasets$nFeature_RNA, prob = 0.01)
merged.datasets <- subset(merged.datasets, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 &
                            percent.mt < 5 & nCount_RNA > countLOW & nCount_RNA < countHIGH)

# Violin plots of all samples
VlnPlot(merged.datasets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Add samples to metadata
sample <- names(merged.datasets@active.ident)
head(sample)
tail(sample)
sample_detect <- ifelse(str_detect(sample, "CRPC1"), "CRPC1", ifelse(str_detect(sample, "CRPC2"), "CRPC2", ifelse(str_detect(sample, "CRPC3"), "CRPC3", ifelse(str_detect(sample, "CRPC4"), "CRPC4", ifelse(str_detect(sample, "CRPC5"), "CRPC5", ifelse(str_detect(sample, "CRPC6"), "CRPC6", ifelse(str_detect(sample, "HEALTH1"), "HEALTH1", ifelse(str_detect(sample, "HEALTH2"), "HEALTH2", ifelse(str_detect(sample, "HEALTH3"), "HEALTH3", ifelse(str_detect(sample, "HEALTH4"), "HEALTH4", ifelse(str_detect(sample, "HEALTH5"), "HEALTH5", ifelse(str_detect(sample, "HEALTH6"), "HEALTH6", ifelse(str_detect(sample, "PCA"), "PCA", NA)))))))))))))
merged.datasets@meta.data$sample <- sample_detect
Idents(object = merged.datasets) <- "sample"
levels(merged.datasets@active.ident) <- unique(sample_detect)

# Normalize, standardize variance, and scale data
options(future.globals.maxSize = 8 * 1024^3)

merged.datasets <- merged.datasets %>%
  SCTransform()
DefaultAssay(merged.datasets) <- "SCT"

#Replace the number of each sample with the group it belongs to:
merged.datasets@meta.data$sample <- recode(merged.datasets@meta.data$sample,
                                           "CRPC1" = "CRPC",
                                           "CRPC2" = "CRPC",
                                           "CRPC3" = "CRPC",
                                           "CRPC4" = "CRPC", 
                                           "CRPC5" = "CRPC",
                                           "CRPC6" = "CRPC",
                                           "HEALTH1" = "HEALTH",
                                           "HEALTH2" = "HEALTH",
                                           "HEALTH3" = "HEALTH",
                                           "HEALTH4" = "HEALTH",
                                           "HEALTH5" = "HEALTH",
                                           "HEALTH6" = "HEALTH",
                                           "PCA" = "PCA")

# Principal Component Analysis
merged.datasets <- RunPCA(merged.datasets, assay= "SCT", npcs = 50)

# Correct batch effects with Harmony
merged.datasets <- RunHarmony(merged.datasets, group.by.vars = "sample")

# Calculate the nearest neighbors matrix
merged.datasets <- FindNeighbors(merged.datasets, reduction = "harmony")

# Cluster
merged.datasets <- FindClusters(merged.datasets, resolution = 1.2)

# Perform UMAP
merged.datasets <- RunUMAP(merged.datasets, reduction = "harmony", dims = 1:40)

# Cluster plot and sample plot
Idents(object = merged.datasets) <- "seurat_clusters"
p1 <- DimPlot(merged.datasets, reduction = "umap", label = TRUE)
p1
Idents(object = merged.datasets) <- "sample"
p2 <- DimPlot(merged.datasets, reduction = "umap", label = TRUE)
p2
p1 + p2

# Save Seurat object
saveRDS(merged.datasets, file = "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat.rds")

# Load Seurat object
#merged.datasets <- readRDS("/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat.rds")

#Violin plot and correlation per group
VlnPlot(merged.datasets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

merged.datasets@meta.data$SCT_snn_res.1.2 <- NULL

# Graphs with top 10 most variable genes
top10 <- head(VariableFeatures(merged.datasets), 10)
top10
FeaturePlot(merged.datasets, features = top10)
RidgePlot(merged.datasets, features = top10)
VlnPlot(merged.datasets, features = top10)
DimPlot(merged.datasets, reduction = "pca")

# Determine the number of PCs
ElbowPlot(merged.datasets)
variance_explained <- merged.datasets[["pca"]]@stdev^2 / sum(merged.datasets[["pca"]]@stdev^2)
cumulative_variance <- cumsum(variance_explained)
plot(cumulative_variance, type = 'b', xlab = 'Number of PCs', ylab = 'Cumulative Variance Explained')

# Perform PCA on the most variable genes
merged.datasets <- RunPCA(merged.datasets, npcs = 30, verbose = FALSE)
merged.datasets <- RunPCA(merged.datasets, features = VariableFeatures(merged.datasets))
print(merged.datasets[["pca"]], dims = 1:10, nfeatures = 5)

#Visualization of PCA
VizDimLoadings(merged.datasets, dims = 1:2, reduction = "pca")
DimPlot(merged.datasets, reduction = "pca")

# Save Seurat object
saveRDS(merged.datasets, file = "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat_2.rds")

# Load Seurat object
#merged.datasets <- readRDS("/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat_2.rds")


# Number of current clusters
Idents(merged.datasets) <- "seurat_clusters"
identifiers <- Idents(merged.datasets)
print(identifiers)

unique_clusters <- unique(identifiers)
print(unique_clusters)

num_clusters <- length(unique_clusters)
print(num_clusters)

table(Idents(merged.datasets))

merged.datasets <- PrepSCTFindMarkers(merged.datasets)
merged.datasets.markers <- FindAllMarkers(merged.datasets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(merged.datasets, file = "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat_3.rds")
write_xlsx(merged.datasets.markers, "/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/FindAllMarkers.xlsx")
#merged.datasets.markers <- read_excel("/home/csic/epi/sbb/Silvia/Analysis/Scripts/FindAllMarkers.xlsx")