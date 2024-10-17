# Load libraries
library(Seurat)
library(Matrix)


# FOR HEALTH SAMPLES:
# Define a function to process health samples
process_seurat <-
  function(file_path) {
  # Read a txt.gz file and creates a count matrix
  count_matrix <- read.table(gzfile(file_path), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Extract row names from the first column
  rownames_temp <- count_matrix[, 1]
  rownames_temp <- make.unique(rownames_temp) # Ensure row names are unique
  # Remove the first colum which contains the row names
  count_matrix <- count_matrix[, -1]
  # Assign the row names to the count matrix
  rownames(count_matrix) <- rownames_temp
  # Convert the count matrix to a standard matrix format
  count_matrix <- as.matrix(count_matrix)
  # Replace underscores in column names with hyphens
  colnames(count_matrix) <- gsub("_", "-", colnames(count_matrix))
  # Extract the sample name from the file name
  sample_name <- gsub("\\.txt(\\.gz)?$", "", basename(file_path))
  # Rename the columns of the count matrix to include the sample name
  colnames(count_matrix) <- paste0(sample_name, "_", seq_len(ncol(count_matrix)))
  # Convert the count matrix to a sparse count matrix
  count_matrix_sparse <- as(count_matrix, "CsparseMatrix")
  # Create a Seurat object from the sparse count matrix
  seurat_object <- CreateSeuratObject(counts = count_matrix_sparse)
  # Normalize, find the variable features, scale the data and PCA
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object)
  # Return the processed Seurat object
  return(seurat_object)
}

# Define the input directory containing the files
input_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/ORIGINAL/HEALTHY"
# List of all files in the input directory that match the pattern *txt.gz
file_list <- list.files(path = input_dir, pattern = "*txt.gz", full.names = TRUE)

# Apply the process_seurat function to each file in file_list
seurat_objects <- lapply(file_list, process_seurat)

# Assign names to each Seurat object based on the file names
names(seurat_objects) <- sapply(file_list,function(x) gsub("\\.txt(\\.gz)?$", "", basename(x)))

# Define the output directory                    
output_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/SAMPLES/"
# Save each Seurat object as an RDS file in the output directory
for (i in seq_along(seurat_objects)) {
  saveRDS(seurat_objects[[i]], file = file.path(output_dir, paste0(names(seurat_objects)[i], ".rds")))
}

#FOR PCA SAMPLES:
# Define a function to process PCA samples
process_seurat <- function(file_path) {
  # Read a gzipped text file into a count matrix
  count_matrix <- read.table(gzfile(file_path), header = TRUE, row.names = 1, sep = "\t")
  # Convert the count matrix to a standard matrix format
  count_matrix <- as.matrix(count_matrix)
  # Replace underscores in column names with hyphens
  colnames(count_matrix) <- gsub("_", "-", colnames(count_matrix))
  # Extract the sample name from the file name
  sample_name <- gsub("\\.txt(\\.gz)?$", "", basename(file_path))
  # Rename the columns of the count matrix to include the sample name
  colnames(count_matrix) <- paste0(sample_name, "_", seq_len(ncol(count_matrix)))
  # Convert the count matrix to a sparse matrix format
  count_matrix_sparse <- as(count_matrix, "CsparseMatrix")
  # Create a Seurat object from the sparse count matrix
  seurat_object <- CreateSeuratObject(counts = count_matrix_sparse)
  # Normalize, FindVariableFeatures, scale data and Run PCA
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object)
  # Return the processed Seurat object
  return(seurat_object)
}

# Define the tar file and the extraction directory                    
tar_file <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/ORIGINAL/PCA.tar"
extract_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/ORIGINAL/"

# Extract the contents of the tar file into the specified directory                               
untar (tar_file, exdir = extract_dir)

# List all gzipped files in the extraction directory                               
file_list <- list.files(path = extract_dir, pattern = "*.txt.gz", full.names = TRUE)

# Apply the process_seurat function to each file in file_list                                
seurat_objects <-lapply(file_list, process_seurat)
# Assign names to each Seurat object based on the file names
names(seurat_objects) <- basename(file_list)

# Define the output directory                                
output_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/SAMPLES/"
# Save each Seurat object as an RDS file in the output directory                                
for (i in seq_along(seurat_objects)) {
  saveRDS(seurat_objects[[i]], file = file.path(output_dir, paste0(names(seurat_objects)[i], ".rds")))
}

                                
#FOR CRPC SAMPLES:
# Define a function to preprocess CRPC samples                             
preprocess <- function(file_path) {
  # Read the existing Seurat object from an RDS file
  seurat_object <- readRDS(file_path)
  # Normalize, find variable features, scale the data and run PCA
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object)
  # Return the processed Seurat object
  return(seurat_object)
}

# Define the input directory where the RDS files are located                                
input_dir <- "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/pre-Seurat/"
# List all RDS files in the input directory                                
file_list <- list.files(path = input_dir, pattern = "*rds", full.names = TRUE)

# Apply the preprocess function to each file in file_list                                
seurat_objects <- lapply(file_list, preprocess)
# Assign names to each Seurat object based on the file names                                
names(seurat_objects) <- basename(file_list)

# Define the output directory where the processed CRPC samples will be saved                                
output_dir <- "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/Seurat"

# Save each processed Seurat object as an RDS file in the output directory                                
for (i in seq_along(seurat_objects)) {
  saveRDS(seurat_objects[[i]], file = file.path(output_dir, paste0(names(seurat_objects)[i], ".rds")))
}
