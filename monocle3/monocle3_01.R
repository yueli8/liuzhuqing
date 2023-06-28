setwd("~/liuzhuqing")
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#https://cole-trapnell-lab.github.io/monocle3/docs/installation/
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="hms_cluster_id_test.rds")
#seurat <- readRDS(file="NKT.rds")
#seurat <- readRDS(file="Treg.rds")
#seurat <- readRDS(file="CD8_T.rds")
#seurat <- readRDS(file="T_helper.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)#pca analysis

#umap,tSNE降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('cds.umap')
p1
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p2
cds <- reduce_dimension(cds, reduction_method = "tSNE")
p3 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="tech")
p3
p4 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype") 
p4

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

colData(cds)$assigned_cell_type <- as.character(partitions(cds))

colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="Mesenchymal",
                                                 "1"="T_helper",
                                                 "2"="B",
                                                 "3"="Mesenchymal",
                                                 "4"="Mesenchymal",
                                                 "5"="Endothelial",
                                                 "6"="CD8_T",
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


cds <- learn_graph(cds)
p1= plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p1
p2 = plot_cells(cds,color_cells_by = "celltype",label_groups_by_cluster = FALSE,label_cell_groups = FALSE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p2


p3 = plot_cells(cds,color_cells_by = "tech",label_groups_by_cluster = FALSE,label_cell_groups = FALSE, 
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p3

p5 = plot_cells(cds,color_cells_by = "assigned_cell_type",label_cell_groups = FALSE, label_groups_by_cluster = TRUE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 6)

p5

##细胞按拟时排序
cds <- order_cells(cds) #dragging a rectangle,then choose,done
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 4)

#NKT cells
p + geom_vline(xintercept = seq(-9,-8,0.25)) + geom_hline(yintercept = seq(-1,0,0.25))
embed <- data.frame(Embeddings(seurat, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -8.25 & UMAP_1 < -8 & UMAP_2 > -0.25 & UMAP_2 < 0)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE, graph_label_size = 6)

#treg cells
p+geom_vline(xintercept = seq(-3,-2,0.25)) + geom_hline(yintercept = seq(-5,-4,0.25))
embed <- data.frame(Embeddings(seurat, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -2.75 & UMAP_1 < -2.5 & UMAP_2 > -5 & UMAP_2 < -4.75)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE, graph_label_size = 5)


##寻找拟时轨迹差异基因
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
#long time
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)#take long time

write.csv(Track_genes,"Trajectory_genes_Treg.csv",row.names = TRUE)
#挑选top10画图展示
#Track_genes_sig <- Track_genes %>% top_n(n=20, morans_I) %>%
 # pull(gene_short_name) %>% as.character()
#基因表达趋势图
#plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="tech", 
#                         min_expr=0.5, ncol = 2)

#plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="celltype", 
#                         min_expr=0.5, ncol = 2)

#plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="seurat_clusters", 
#                         min_expr=0.5, ncol = 2)

#plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="assigned_cell_type", 
#                         min_expr=0.5, ncol = 2)#plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="assigned_cell_type", 
#                         min_expr=0.5, ncol = 2)
#FeaturePlot图
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)

p$facet$params$ncol <- 5

#13 genes    morans_I>0.5 compared with 167_immune_gene
Track<-c("CXCR4","CD3E","CD163","CD3D","MRC1","HLA-DRB1","C5AR1","FCGR2A","IL7R","IL32","CD69","CD14","LYZ")

#NKT   8 genes    morans_I>0.1 compared with 167_immune_gene
#Track<-c("CD8B","TNFRSF18","IL32","FCGR3A","CD3D","CD3G","GZMH","CCL5")

#Treg  5 genes    morans_I>0.1 compared with 167_immune_gene
#Track<-c("BIRC3","CD69","HLA-A","HLA-DRB1","MAL")

#CD8  10 genes    morans_I>0.1 compared with 167_immune_gene
#Track<-c("CCL5","CD3E","CD69","CXCR4","GZMK","ICAM1","IL2RA","TNFRSF18","TNFRSF4","XCL1")

#Thelper 1 genes    morans_I>0.1 compared with 167_immune_gene
#Track<-c("CXCR4")


plot_genes_in_pseudotime(cds[Track,], color_cells_by="tech", 
                         min_expr=0.5, ncol = 2)

plot_genes_in_pseudotime(cds[Track,], color_cells_by="celltype", 
                         min_expr=0.5, ncol = 2)

plot_genes_in_pseudotime(cds[Track,], color_cells_by="seurat_clusters", 
                         min_expr=0.5, ncol = 2)

plot_genes_in_pseudotime(cds[Track,], color_cells_by="assigned_cell_type", 
                         min_expr=0.5, ncol = 2)




#FeaturePlot图
Track<-c("CXCR4","CD3D")

plot_cells(cds, genes=Track, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)

NKT<-readRDS("NKT.rds")

markers.to.plot<-c("CXCR4","CD3E","CD163","CD3D","MRC1","HLA-DRB1","C5AR1","FCGR2A","IL7R","CD69","CD14")


plot_cells(NKT, genes=Track, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)


##寻找共表达模块(work)
gene_short_name = rownames(data)
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)#take long time
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
#gene_module <- find_gene_modules(cds[genelist,], resolution=1e-1, cores = 6)
write.csv(gene_module,"Genes_Module.csv",row.names=F)

cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$predicted.id)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group = NULL)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p<-pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
p




#heatmap
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
library(pheatmap)
#take long time
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

write.csv(modulated_genes,"modulated_genes.csv",row.names = TRUE)

genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.52))
genes

pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 6 groups

pheatmap(pt.matrix)

pheatmap(pt.matrix, cellheight=10, cluster_cols=T,cluster_rows = T, 
         color=colorRampPalette(c("green", "black", "red"))(100),fontsize=7)

head(cds@metadata)
colnames(colData(cds))


write.csv(assays(cds),"assays")
a1<-read.table("assays",sep=",",row.names = 1)
select_gene<-c("CCL5","GZMA","CCL5","GZMA","NKG7","CCL4","CST3","CD68","TYROBP","FCER1G","APOC1","HBB")
a2<-a1[select_gene,]
a3<-t(a2)
a4<-a3[,1]

