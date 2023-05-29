setwd("~/liuzhuqing")
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
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p

#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="Mesenchymal",
                                                 "1"="NKT",
                                                 "2"="B",
                                                 "3"="Mesenchymal",
                                                 "4"="Mesenchymal",
                                                 "5"="Endothelial",
                                                 "6"="NKT",
                                                 "7"="Mesenchymal",
                                                 "8"="Th17",
                                                 "9"="Macrophage",
                                                 "10"="Epithelial",
                                                 "11"="Dendritic",
                                                 "12"="Mesenchymal",
                                                 "13"="NKT",
                                                 "14"="Endothelial",
                                                 "15"="Basal",
                                                 "16"="Monocyte",
                                                 "17"="Treg",
                                                 "18"="Treg")


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
plotCytoTRACE(result01, phenotype = phe,gene = "CXCR4",outputDir = "cxcr4")
plotCytoTRACE(result01, phenotype = phe,gene = "CD3E",outputDir = "cd3e")
plotCytoTRACE(result01, phenotype = phe,gene = "CD163",outputDir = "cd163")
plotCytoTRACE(result01, phenotype = phe,gene = "CD3D",outputDir = "cd3d")
plotCytoTRACE(result01, phenotype = phe,gene = "MRC1",outputDir = "mrc1")
plotCytoTRACE(result01, phenotype = phe,gene = "HLA-DRB1",outputDir = "hla-drb1")
plotCytoTRACE(result01, phenotype = phe,gene = "CSAR1",outputDir = "csar1")
plotCytoTRACE(result01, phenotype = phe,gene = "FCGR2A",outputDir = "fcgr2a")
plotCytoTRACE(result01, phenotype = phe,gene = "IL7R",outputDir = "il7r")
plotCytoTRACE(result01, phenotype = phe,gene = "CD69",outputDir = "cd69")
plotCytoTRACE(result01, phenotype = phe,gene = "CD14",outputDir = "cd14")
plotCytoTRACE(result01, phenotype = phe,gene = "BIRC3",outputDir = "birc3")
plotCytoTRACE(result01, phenotype = phe,gene = "CCL5",outputDir = "ccl5")
plotCytoTRACE(result01, phenotype = phe,gene = "CD3G",outputDir = "cd3g")
plotCytoTRACE(result01, phenotype = phe,gene = "CD8B",outputDir = "cd8b")
plotCytoTRACE(result01, phenotype = phe,gene = "FCGR3A",outputDir = "fcgr3a")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMK",outputDir = "gzmk")
plotCytoTRACE(result01, phenotype = phe,gene = "GZMH",outputDir = "gzmh")
plotCytoTRACE(result01, phenotype = phe,gene = "HLA-A",outputDir = "hla-a")
plotCytoTRACE(result01, phenotype = phe,gene = "ICAM1",outputDir = "icam1")
plotCytoTRACE(result01, phenotype = phe,gene = "IL2RA",outputDir = "il2ra")
plotCytoTRACE(result01, phenotype = phe,gene = "IL32",outputDir = "il32")
plotCytoTRACE(result01, phenotype = phe,gene = "LYZ",outputDir = "lyz")
plotCytoTRACE(result01, phenotype = phe,gene = "MAL",outputDir = "mal")
plotCytoTRACE(result01, phenotype = phe,gene = "TNFRSF18",outputDir = "tnfrsf18")
plotCytoTRACE(result01, phenotype = phe,gene = "TNFRSF4",outputDir = "tnfrsf4")
plotCytoTRACE(result01, phenotype = phe,gene = "XCL1",outputDir = "xcl1")


