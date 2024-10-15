## CREATE SEURAT OBJECT POR CRPC:

library(dplyr)
library(Seurat)
library(tidyverse)

#File 1:
#CRPC_1 <- read_csv("/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428952_1778_JZ_HMP_04_IGO_10726_2_dense.csv")

# Eliminar la primera columna y la de clusters, pero meter la info de clusters en una nueva variable
#cluster_info <- CRPC_1$CLUSTER
#count_matrix <- CRPC_1 %>% select(-c(...1, CLUSTER))
#count_matrix <- as.data.frame(count_matrix)

# Convertir la matriz en un df con nombres de filas (números de las células)
#Aquí cambia CRPC1 para cada paciente
#rownames(count_matrix) <- paste0("CRPC1_", seq_len(nrow(count_matrix)))

# Crear el objeto Seurat
#seurat_CRPC1 <- CreateSeuratObject(counts = t(count_matrix), project = "scRNAseq")

# Añadir los clústeres como metadatos
#seurat_CRPC1$Cluster <- cluster_info

#Guardar Seurat:
#saveRDS(seurat_CRPC1, file = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/CRPC_1.rds")


#REPETIMOS PARA RESTO DE MATRICES:
#CRPC_2 <- read_csv("/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428953_1779_HMP05_IGO_10726_3_dense.csv")

#cluster_info2 <- CRPC_2$CLUSTER
#count_matrix2 <- CRPC_2 %>% select(-c(...1, CLUSTER))
#count_matrix2 <- as.data.frame(count_matrix2)

#rownames(count_matrix2) <- paste0("CRPC2_", seq_len(nrow(count_matrix2)))

#seurat_CRPC2 <- CreateSeuratObject(counts = t(count_matrix2), project = "scRNAseq")

#seurat_CRPC2$Cluster <- cluster_info2

#saveRDS(seurat_CRPC2, file = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/CRPC_2.rds")

#CRPC_3:
#CRPC_3 <- read_csv("/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428954_1845_HMP-08_IGO_10837_19_dense.csv")

#cluster_info3 <- CRPC_3$CLUSTER
#count_matrix3 <- CRPC_3 %>% select(-c(...1, CLUSTER))
#count_matrix3 <- as.data.frame(count_matrix3)

#rownames(count_matrix3) <- paste0("CRPC3_", seq_len(nrow(count_matrix3)))

#seurat_CRPC3 <- CreateSeuratObject(counts = t(count_matrix3), project = "scRNAseq")

#seurat_CRPC3$Cluster <- cluster_info3

#saveRDS(seurat_CRPC3, file = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/CRPC_3.rds")

#CRPC_4:
#CRPC_4 <- read_csv("/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428955_1968_HMP11_1_IGO_11247_3_dense.csv")

#cluster_info4 <- CRPC_4$CLUSTER
#count_matrix4 <- CRPC_4 %>% select(-c(...1, CLUSTER))
#count_matrix4 <- as.data.frame(count_matrix4)

#rownames(count_matrix4) <- paste0("CRPC4_", seq_len(nrow(count_matrix4)))

#seurat_CRPC4 <- CreateSeuratObject(counts = t(count_matrix4), project = "scRNAseq")

#seurat_CRPC4$Cluster <- cluster_info4

#saveRDS(seurat_CRPC4, file = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/pre-Seurat/CRPC_4.rds")


#CRPC_:
#CRPC_5 <- read_csv("/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428956_1969_HMP11_2_IGO_11247_4_dense.csv")

#cluster_info5 <- CRPC_5$CLUSTER
#count_matrix5 <- CRPC_5 %>% select(-c(...1, CLUSTER))
#count_matrix5 <- as.data.frame(count_matrix5)

#rownames(count_matrix5) <- paste0("CRPC5_", seq_len(nrow(count_matrix5)))

#seurat_CRPC5 <- CreateSeuratObject(counts = t(count_matrix5), project = "scRNAseq")

#seurat_CRPC5$Cluster <- cluster_info5

#saveRDS(seurat_CRPC5, file = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/CRPC_5.rds")

#CRPC_6:
#CRPC_6 <- read_csv("/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/GSM6428957_1846_JZHP_3_IGO_10837_21_dense.csv")

#cluster_info6 <- CRPC_6$CLUSTER
#count_matrix6 <- CRPC_6 %>% select(-c(...1, CLUSTER))
#count_matrix6 <- as.data.frame(count_matrix6)

#rownames(count_matrix6) <- paste0("CRPC6_", seq_len(nrow(count_matrix6)))

#seurat_CRPC6 <- CreateSeuratObject(counts = t(count_matrix6), project = "scRNAseq")

#seurat_CRPC6$Cluster <- cluster_info6

#saveRDS(seurat_CRPC6, file = "/home/csic/epi/sbb/Silvia/Analysis/Datos/Zaidi et al./Matrix/CRPC/CRPC_6.rds")