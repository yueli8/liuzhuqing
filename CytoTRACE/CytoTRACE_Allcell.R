setwd("~/liuzhuqing/herong")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="hms_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p
#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="Endothelial",  "1"="Myofibroblast", "2"="Inflammatory_Fibroblast",
                                                 "3"="Melanocyte", "4"="Inflammatory_Fibroblast", "5"="Endothelial",
                                                 "6"="T",   "7"="Inflammatory_Fibroblast", "8"="Myofibroblast",  "9"="T",
                                                 "10"="Dendritic",  "11"="Myofibroblast", "12"="Inflammatory_Fibroblast",
                                                 "13"="Melanocyte", "14"="Inflammatory_Fibroblast", "15"="Epithelial",
                                                 "16"="Epithelial",   "17"="Natural_Killer", "18"="Dendritic","19"="Melanocyte",
                                                 "20"="Inflammatory_Fibroblast",  "21"="Endothelial", "22"="Epithelial",
                                                 "23"="Melanocyte", "24"="Mast", "25"="Neutrophil","26"="Epithelial", 
                                                 "27"="Endothelial", "28"="Melanocyte", "29"="Inflammatory_Fibroblast", 
                                                 "30"="Inflammatory_Fibroblast",  "31"="Melanocyte", "32"="Epithelial")
colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)

setwd("~/liuzhuqing/CytoTRACE")
plotCytoTRACE(result01, phenotype = phe)
plotCytoTRACE(result01, phenotype = phe,gene = "C5AR1",outputDir = "c5ar1")
plotCytoTRACE(result01, phenotype = phe,gene = "CCL5",outputDir = "ccl5")
plotCytoTRACE(result01, phenotype = phe,gene = "CD163",outputDir = "cd163")
plotCytoTRACE(result01, phenotype = phe,gene = "CD3D",outputDir = "cd3d")
plotCytoTRACE(result01, phenotype = phe,gene = "CD3E",outputDir = "cd3e")
plotCytoTRACE(result01, phenotype = phe,gene = "CXCR4",outputDir = "cxcr4")
plotCytoTRACE(result01, phenotype = phe,gene = "FCGR3A",outputDir = "fcgr3a")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMA",outputDir = "gzma")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMB",outputDir = "gzmb")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMH",outputDir = "gzmh")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMK",outputDir = "gzmk")
plotCytoTRACE(result01, phenotype = phe,gene = "HLA-DRB1",outputDir = "hla-drb1")
plotCytoTRACE(result01, phenotype = phe,gene = "IL32",outputDir = "il32")
plotCytoTRACE(result01, phenotype = phe,gene = "IL7R",outputDir = "il7r")
plotCytoTRACE(result01, phenotype = phe,gene = "LYZ",outputDir = "lyz")
plotCytoTRACE(result01, phenotype = phe,gene = "MRC1",outputDir = "mrc1")
plotCytoTRACE(result01, phenotype = phe,gene = "PRF1",outputDir = "PRF1")
plotCytoTRACE(result01, phenotype = phe,gene = "SOX10",outputDir = "sox10")