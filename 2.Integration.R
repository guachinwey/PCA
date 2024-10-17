# Load libraries
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

# Define input directory
input_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/SAMPLES/"
# List .rds files in the directory
file_list <- list.files(path = input_dir, pattern = "*.rds", full.names = TRUE)

# Define names for the Seurat objects corresponding to the samples
# Assign names for each Seurat object that wil be loaded later
object_names <- c("CRPC_1", "CRPC_2", "CRPC_3", "CRPC_4", "CRPC_5", "CRPC_6", 
                  "HEALTH1", "HEALTH2", "HEALTH3", "HEALTH4", "HEALTH5", "HEALTH6", "PCA")

# Load Seurat objects and assign names
# Load each .rds file as a Seurat object and assign them to the defined names
seurat_objects <- lapply(file_list, readRDS)
names(seurat_objects) <- object_names

# Rename the active identity levels of each object to match their sample name
# This ensures that the identity levels correspond to the sample name
for (i in 1:13) {
  levels(seurat_objects[[i]]@active.ident) <- object_names[i]
}
# Assign the objects to the global environment for easier access
list2env(seurat_objects, .GlobalEnv)

# Object verification and dimensions check
# Check the numer of objects and their dimensions for verification
length(seurat_objects) # Print number of objects
lapply(seurat_objects, function(x) levels(x@active.ident)) # Print identity levels for each object
lapply(seurat_objects, function(x) dim(x)) # Print dimensions of each object

# Merge all datasets into one Seurat object, keeping the project name
merged.datasets <- merge(seurat_objects[[1]], y = seurat_objects[-1], project = "PCa_Human")
       
# Display the head of column names, row names and activity identities to verify merging
head(colnames(merged.datasets))
head(merged.datasets)
head(rownames(merged.datasets))
table(merged.datasets@active.ident) # Check the distribution of samples in merged objects
ncol(merged.datasets) # Print number of features in the merged dataset

# Pre-processing: Calculate percentage of mitochondrial genes
# Calculate the percentage of mitochondrial gene expression to assess cell quality
merged.datasets[["percent.mt"]] <- PercentageFeatureSet(merged.datasets, pattern = "^MT-")

# Generate violin plot to assess relationships between features 
# Visualize the distribution of features to verify data quality
VlnPlot(merged.datasets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 # Combine scatter plots

# Set thesholds for filtering cells based on low and high counts/features
# Define criteria to remove low-quality cells
minCov <- 1000
countLOW <- ifelse(min(merged.datasets$nCount_RNA) >= minCov, min(merged.datasets$nCount_RNA), quantile(merged.datasets$nCount_RNA, prob = c(0.01)))
countHIGH <- quantile(merged.datasets$nCount_RNA, prob = 0.99)
featureLOW <- quantile(merged.datasets$nFeature_RNA, prob = 0.01)

# Subset the merged dataset to remove low-quality cells
# Apply the filtering criteria to remove cells that do not meet the thresholds
merged.datasets <- subset(merged.datasets, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 &
                            percent.mt < 5 & nCount_RNA > countLOW & nCount_RNA < countHIGH)

# Visualize the filtered features using violin plots again
# Visually check the features after filtering
VlnPlot(merged.datasets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Add samples information to the metadata for further analysis
# Detect and assign sample groups based on identifiers in the sample names
sample <- names(merged.datasets@active.ident)

# Sample detection using regular expressions
# Identify and categorize samples according to their names
sample_detect <- ifelse(str_detect(sample, "CRPC1"), "CRPC1", 
                        ifelse(str_detect(sample, "CRPC2"), "CRPC2", 
                               ifelse(str_detect(sample, "CRPC3"), "CRPC3", 
                                      ifelse(str_detect(sample, "CRPC4"), "CRPC4", 
                                             ifelse(str_detect(sample, "CRPC5"), "CRPC5", 
                                                    ifelse(str_detect(sample, "CRPC6"), "CRPC6", 
                                                           ifelse(str_detect(sample, "HEALTH1"), "HEALTH1", 
                                                                  ifelse(str_detect(sample, "HEALTH2"), "HEALTH2", 
                                                                         ifelse(str_detect(sample, "HEALTH3"), "HEALTH3", 
                                                                                ifelse(str_detect(sample, "HEALTH4"), "HEALTH4", 
                                                                                       ifelse(str_detect(sample, "HEALTH5"), "HEALTH5", 
                                                                                              ifelse(str_detect(sample, "HEALTH6"), "HEALTH6", 
                                                                                                     ifelse(str_detect(sample, "PCA"), "PCA", NA)))))))))))))
merged.datasets@meta.data$sample <- sample_detect # Add detected samples to metadata
Idents(object = merged.datasets) <- "sample" # Set identity based on samples
levels(merged.datasets@active.ident) <- unique(sample_detect) # Update levels

# Normalize the data and perform variance stabilization and scaling
# Normalize using SCTransform to stabilize the variance
options(future.globals.maxSize = 8 * 1024^3) # Increase memory limit if needed
merged.datasets <- merged.datasets %>%
  SCTransform() # Use SCTransform to normalize and variance stabilize
DefaultAssay(merged.datasets) <- "SCT" # Set the default assay to the normalized data

# Replace sample identifiers with broader group labels (CRPC, HEALTH and PCA)
# Simplification of sample identifiers to more general categories
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

# Perform Principal Component Analysis (PCA) on the normalized data
# Apply PCA to reduce the dimensionality of the normalized data
merged.datasets <- RunPCA(merged.datasets, assay= "SCT", npcs = 50)

# Correct for batch effects using Harmony based on the sample groups
# Adjust the data structure to remove unwanted batch effects
merged.datasets <- RunHarmony(merged.datasets, group.by.vars = "sample")

# Calculate the nearest neighbors for clustering
# Prepare the object for clustering analysis
merged.datasets <- FindNeighbors(merged.datasets, reduction = "harmony")

# Perform clustering on the data
merged.datasets <- FindClusters(merged.datasets, resolution = 1.2)
       
# Run UMAP for dimensionality reduction and visualization
# This step reduces the high-dimensional data to two dimension for easier visualization
merged.datasets <- RunUMAP(merged.datasets, reduction = "harmony", dims = 1:40)

# Generate UMAP plots to visualize clusters and sample distribution
# Set the identities to the clusters for visualization
Idents(object = merged.datasets) <- "seurat_clusters" # Set identities to clusters
p1 <- DimPlot(merged.datasets, reduction = "umap", label = TRUE) # UMAP plot for clusters
p1 
Idents(object = merged.datasets) <- "sample" # Set identities to sample groups
p2 <- DimPlot(merged.datasets, reduction = "umap", label = TRUE) # UMAP plot for samples
p2
p1 + p2 # Combine the two UMAP plots

# Save the merged Seurat object for future analysis
saveRDS(merged.datasets, file = "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat.rds")

# Load the Seurat object if needed for further analysis
#merged.datasets <- readRDS("/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat.rds")

# Generate violin plots and scatter plots to assess quality metrics
# Visualize the distribution of features in the dataset to assess data quality
VlnPlot(merged.datasets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Create scatter plots to explore the relationships between features
plot1 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "percent.mt") # Count vs. Mitochondrial percentages
plot2 <- FeatureScatter(merged.datasets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") # Count vs. Number of features
plot1 + plot2 # Combine the scatter plots for visual comparison

# Remove the clustering result from the metadata       
merged.datasets@meta.data$SCT_snn_res.1.2 <- NULL

# Identify and vidualize the top 10 most variable genes in the dataset
top10 <- head(VariableFeatures(merged.datasets), 10) # Get the top 10 variable genes
top10 

# Create various plots for the top variable genes to visualize their expression across clusters
# Visualize the top variable genes using different plotting methods
FeaturePlot(merged.datasets, features = top10) # Feature plot for the top variable genes
RidgePlot(merged.datasets, features = top10) # Ridge plot for the top variable genes
VlnPlot(merged.datasets, features = top10) # Violin plot for the top variable genes
DimPlot(merged.datasets, reduction = "pca") # PCA plot to visualize samples in PCA space

# Determine the optimal number of principa components (PCs) 
ElbowPlot(merged.datasets) # Plot to visualize the variance explained by each PC
# Calculate the variance explaines by each PC
variance_explained <- merged.datasets[["pca"]]@stdev^2 / sum(merged.datasets[["pca"]]@stdev^2)
# Calculate cumulative variance for the PCs
cumulative_variance <- cumsum(variance_explained)
# Plot cumulative variance to assess how many PCs to keep
plot(cumulative_variance, type = 'b', xlab = 'Number of PCs', ylab = 'Cumulative Variance Explained')

# Perform PCA specifically on the most variable genes identified earlier
merged.datasets <- RunPCA(merged.datasets, npcs = 15, verbose = FALSE)
# Run PCA again focusing on the most variable genes
merged.datasets <- RunPCA(merged.datasets, features = VariableFeatures(merged.datasets))
# Print the PCA results for the top 10 dimensions and their associated features
print(merged.datasets[["pca"]], dims = 1:10, nfeatures = 5)

# Visualization of PCA loading for the firts two dimensions to understand gene contributions
# Visualize how genes contibute to the first two principal components
VizDimLoadings(merged.datasets, dims = 1:2, reduction = "pca")
# Plot PCA results to visualize the distribution of samples in PCA space
DimPlot(merged.datasets, reduction = "pca")

# Save the Seurat object
saveRDS(merged.datasets, file = "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat_2.rds")

# Load Seurat object if needed
#merged.datasets <- readRDS("/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat_2.rds")

# Check the number of clusters identified in the dataset
Idents(merged.datasets) <- "seurat_clusters" # Set the identity to the cluster IDs
identifiers <- Idents(merged.datasets) # Get the cluster identifiers
print(identifiers) # Print clusters identifiers for review

# Get unique cluster identifiers to determine how many clusters exist       
unique_clusters <- unique(identifiers) 
print(unique_clusters)

# Count the number of clusters detected
num_clusters <- length(unique_clusters)
print(num_clusters)

# Generate a frequency table to see how many cells belong to each cluster
table(Idents(merged.datasets))

# Prepare the dataset for marker detection using SCT method      
# Set up the object for detecting marker genes between clusters
merged.datasets <- PrepSCTFindMarkers(merged.datasets) # Prepare for finding markers
# Find markets for all clusters compared to all other clusters
# Identify positive markers with specific thresholds for inclusion
merged.datasets.markers <- FindAllMarkers(merged.datasets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # Filtering to only include positive markerks with minimum percentage and log-fold change thresholds 

# Save the Seurat object
saveRDS(merged.datasets, file = "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/INTEGRATION/harmony_seurat_3.rds")

# Write the markers found into an Excel file for external analysis and reporting
write_xlsx(merged.datasets.markers, "/home/csic/epi/sbb/Silvia/Analysis/NEWS/RESULTS/FindAllMarkers.xlsx")

# Uncomment the following line to read markers from the saves Excel file if needed
#merged.datasets.markers <- read_excel("/home/csic/epi/sbb/Silvia/Analysis/Scripts/FindAllMarkers.xlsx")
