library(Seurat)
library(Matrix)


#FOR HEALTH:
process_seurat <-
  function(file_path) {
  count_matrix <- read.table(gzfile(file_path), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames_temp <- count_matrix[, 1]
  rownames_temp <- make.unique(rownames_temp)
  count_matrix <- count_matrix[, -1]
  rownames(count_matrix) <- rownames_temp
  count_matrix <- as.matrix(count_matrix)
  colnames(count_matrix) <- gsub("_", "-", colnames(count_matrix))
  sample_name <- gsub("\\.txt(\\.gz)?$", "", basename(file_path))
  colnames(count_matrix) <- paste0(sample_name, "_", seq_len(ncol(count_matrix)))
  count_matrix_sparse <- as(count_matrix, "CsparseMatrix")
  seurat_object <- CreateSeuratObject(counts = count_matrix_sparse)
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object)
  return(seurat_object)
}

input_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/ORIGINAL/HEALTHY"
file_list <- list.files(path = input_dir, pattern = "*txt.gz", full.names = TRUE)
seurat_objects <- lapply(file_list, process_seurat)
names(seurat_objects) <- sapply(file_list,function(x) gsub("\\.txt(\\.gz)?$", "", basename(x)))

output_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/SAMPLES/"
for (i in seq_along(seurat_objects)) {
  saveRDS(seurat_objects[[i]], file = file.path(output_dir, paste0(names(seurat_objects)[i], ".rds")))
}

#For PCA:
process_seurat <- function(file_path) {
  count_matrix <- read.table(gzfile(file_path), header = TRUE, row.names = 1, sep = "\t")
  count_matrix <- as.matrix(count_matrix)
  colnames(count_matrix) <- gsub("_", "-", colnames(count_matrix))
  sample_name <- gsub("\\.txt(\\.gz)?$", "", basename(file_path))
  colnames(count_matrix) <- paste0(sample_name, "_", seq_len(ncol(count_matrix)))
  count_matrix_sparse <- as(count_matrix, "CsparseMatrix")
  seurat_object <- CreateSeuratObject(counts = count_matrix_sparse)
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object)
  return(seurat_object)
}

tar_file <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/ORIGINAL/PCA.tar"
extract_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/ORIGINAL/"

untar (tar_file, exdir = extract_dir)

file_list <- list.files(path = extract_dir, pattern = "*.txt.gz", full.names = TRUE)
seurat_objects <-lapply(file_list, process_seurat)
names(seurat_objects) <- basename(file_list)

output_dir <- "/home/csic/epi/sbb/Silvia/Analysis/NEWS/DATA/SEURAT/SAMPLES/"
for (i in seq_along(seurat_objects)) {
  saveRDS(seurat_objects[[i]], file = file.path(output_dir, paste0(names(seurat_objects)[i], ".rds")))
}

#FOR CRPC:
#preprocess <- function(file_path) {
#  seurat_object <- readRDS(file_path)
#  seurat_object <- NormalizeData(seurat_object)
#  seurat_object <- FindVariableFeatures(seurat_object)
#  seurat_object <- ScaleData(seurat_object)
#  seurat_object <- RunPCA(seurat_object)
#  return(seurat_object)
#}
#
#input_dir <- "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/pre-Seurat/"
#file_list <- list.files(path = input_dir, pattern = "*rds", full.names = TRUE)
#seurat_objects <- lapply(file_list, preprocess)
#names(seurat_objects) <- basename(file_list)
#
#output_dir <- "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/Seurat"
#for (i in seq_along(seurat_objects)) {
#  saveRDS(seurat_objects[[i]], file = file.path(output_dir, paste0(names(seurat_objects)[i], ".rds")))
#}