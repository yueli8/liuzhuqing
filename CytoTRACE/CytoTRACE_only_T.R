setwd("~/liuzhuqing/herong")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="t.rds")
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
                                                 "0"="T_helper",  "1"="CD4_T", "2"="Treg",
                                                 "3"="CD8_T", "4"="Naive_T", "5"="CD8_T",
                                                 "6"="MAIT",   "7"="Naive_T", "8"="CD8_T",  "9"="NKT",
                                                 "10"="Treg",  "11"="CD8_T", "12"="CD8_T")
colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)

setwd("~/liuzhuqing/CytoTRACE_T")
plotCytoTRACE(result01, phenotype = phe)
plotCytoTRACE(result01, phenotype = phe,gene = "CCL5",outputDir = "ccl5")
plotCytoTRACE(result01, phenotype = phe,gene = "CD8A",outputDir = "cd8a")
plotCytoTRACE(result01, phenotype = phe,gene = "FOXP3",outputDir = "foxp3")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMA",outputDir = "gzma")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMB",outputDir = "gzmb")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMH",outputDir = "gzmh")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMK",outputDir = "gzmk")
plotCytoTRACE(result01, phenotype = phe,gene = "TNFRSF18",outputDir = "tnfrsf18")
plotCytoTRACE(result01, phenotype = phe,gene = "XCL1",outputDir = "xcl1")