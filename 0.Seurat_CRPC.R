# Load libraries
library(dplyr)
library(Seurat)
library(tidyverse)

# Define a function to create a Seurat object from a given file path
create_seurat_object <- function(file_path, sample_name) {
  # Read the CSV file
  data <- read_csv(file_path)
  # Extract cluster information
  cluster_info <- data$CLUSTER
  # Remove the first column and the CLUSTER column to create the count matrix
  count_matrix <- data %>% select(-c(...1, CLUSTER))
  # Convert the count matrix to a data frame
  count_matrix <- as.data.frame(count_matrix)
  # Assign unique row names based on the sample name
  rownames(count_matrix) <- paste0(sample_name, "_", seq_len(nrow(count_matrix)))
  # Create the Seurat object, transposing the count matrix
  seurat_object <- CreateSeuratObject(counts = t(count_matrix), project = "scRNAseq")
  # Add the cluster information as metadata
  seurat_object$Cluster <- cluster_info
  # Return the Seurat object
  return(seurat_object)
}

# Define a list of file paths and corresponding sample names
file_paths <- list(
  CRPC_1 = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428952_1778_JZ_HMP_04_IGO_10726_2_dense.csv",
  CRPC_2 = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428953_1779_HMP05_IGO_10726_3_dense.csv",
  CRPC_3 = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428954_1845_HMP-08_IGO_10837_19_dense.csv",
  CRPC_4 = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428955_1968_HMP11_1_IGO_11247_3_dense.csv",
  CRPC_5 = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428956_1969_HMP11_2_IGO_11247_4_dense.csv",
  CRPC_6 = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428957_1846_JZHP_3_IGO_10837_21_dense.csv"
)

# Create an empty list to store the Seurat objects
seurat_objects <- list()

# Loop through the file paths and create Seurat objects
for (sample in names(file_paths)) {
  seurat_objects[[sample]] <- create_seurat_object(file_paths[[sample]], sample)
  
  # Save each Seurat object as an RDS file
  saveRDS(seurat_objects[[sample]], 
          file = file.path("/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/pre-Seurat/", paste0(sample, ".rds")))
}
