library(Seurat)
library(scRepertoire)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(dplyr)
library(Matrix)
library(muscat)
library(reshape2)
library(celldex)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(DEsingle)
library(stringr)
library(sva)
library(readxl)
library(DESeq2)
library(DESeq)
library(pamr)
library(ggpubr)
library(ggraph)
library(gplots)
library(pca3d)
library(rgl)
library(scatterplot3d)
library(FactoMineR)
library(ggfortify)
library(useful)
library(tidyverse)
library(kableExtra)
library(xfun)
library(psych)
library(limma)
library(calibrate)
library(pheatmap)
library(ggraph)
library(circlize)
library(scales)
#devtools::install_github("ncborcherding/scRepertoire@dev")
#manual: https://ncborcherding.github.io/vignettes/vignette.html
#not install from bioconductor,it is from github. 

#setwd("~/liuzhuqing/tcr")
#tumor and para samples
#Tumor <- read.csv("t_total.csv")
#Para <- read.csv("p_total.csv")
#contig_list <- list(Tumor,Para)
#combined <- combineTCR(contig_list, 
#                       samples = c("Tumor", "Para"), 
#                       ID = c("T", "P"), cells ="T-AB")

setwd("~/liuzhuqing/tcr")
t1<- read.csv("t1_filtered_contig_annotations.csv")
t2<- read.csv("t2_filtered_contig_annotations.csv")  
t3<- read.csv("t3_filtered_contig_annotations.csv")
p1<- read.csv("p1_filtered_contig_annotations.csv")
p2<- read.csv("p2_filtered_contig_annotations.csv")
p3<- read.csv("p3_filtered_contig_annotations.csv")

contig_list <- list(t1,t2,t3,p1,p2,p3)
combined <- combineTCR(contig_list, 
                       samples = c("t1", "t2","t3","p1", "p2", "p3"), 
                       ID = c( "t","t","t","p", "p", "p"), cells ="T-AB")

head(contig_list[[1]])

quantContig(combined, cloneCall="gene+nt", scale = T)
quantContig(combined, cloneCall="gene+nt", chain = "TRA")
quantContig(combined, cloneCall="gene+nt", chain = "TRB")
quantContig(combined, cloneCall="gene+nt", chain = "TRA", scale = TRUE)
quantContig(combined, cloneCall="gene+nt", scale = T, chain = "both")
quantContig_output <- quantContig(combined, cloneCall="gene+nt", 
                                  scale = T, exportTable = T)
quantContig_output

abundanceContig(combined, cloneCall = "gene", scale = F)
abundanceContig(combined, cloneCall = "gene", exportTable = T)
lengthContig(combined, cloneCall="aa", chain = "both") 
lengthContig(combined, cloneCall="nt", chain = "TRA") 

# t1 should be t1_t,not only t1
compareClonotypes(combined, 
                  numbers = 5, 
                  samples = c("t1_t","p1_p"), 
                  cloneCall="aa", 
                  graph = "alluvial")

#export table, no graph
compareClonotypes(combined, 
                numbers = 5, 
                samples = c("t1_t","p1_p"), 
                cloneCall="aa", 
                graph = "alluvial",exportTable = TRUE)

compareClonotypes(combined, 
                  numbers = 5, 
                  samples = c("t2_t","p2_p"), 
                  cloneCall="aa", 
                  graph = "alluvial")

#export table, no graph
compareClonotypes(combined, 
                  numbers = 5, 
                  samples = c("t2_t","p2_p"), 
                  cloneCall="aa", 
                  graph = "alluvial",exportTable = TRUE)

compareClonotypes(combined, 
                  numbers = 5, 
                  samples = c("t3_t","p3_p"), 
                  cloneCall="aa", 
                  graph = "alluvial")

#export table, no graph
compareClonotypes(combined, 
                  numbers = 5, 
                  samples = c("t3_t","p3_p"), 
                  cloneCall="aa", 
                  graph = "alluvial",exportTable = TRUE)

#only one common aa, no good resullt
compareClonotypes(combined, 
                  numbers = 4,      
                  samples = c("t1_t","p1_p","t2_t","p2_p","t3_t","p3_p"), 
                  cloneCall="aa", 
                  graph = "alluvial")

compareClonotypes(combined, 
                  numbers = 4,      
                  samples = c("t1_t","p1_p","t2_t","p2_p","t3_t","p3_p"), 
                  cloneCall="aa", 
                  graph = "alluvial",exportTable = TRUE)

#visualize gene usage
vizGenes(combined, gene = "V", 
         chain = "TRA", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)
vizGenes(combined, gene = "V", 
         chain = "TRB", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)
vizGenes(combined[c(1,2,3)], 
         gene = "V", 
         chain = "TRB", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")
vizGenes(combined[c(4,5,6)], 
         gene = "V", 
         chain = "TRB", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")

vizGenes(combined[c(1,2,3)], 
         gene = "V", 
         chain = "TRA", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")
vizGenes(combined[c(4,5,6)], 
         gene = "V", 
         chain = "TRA", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")

clonalHomeostasis(combined, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))

clonalHomeostasis(combined, cloneCall = "aa")
clonalProportion(combined, cloneCall = "gene",
                 split = c(10, 100, 1000, 10000, 30000, 1e+05)) 
clonalProportion(combined, cloneCall = "nt") 
clonalOverlap(combined, cloneCall = "gene+nt", 
              method = "morisita")
clonesizeDistribution(combined, cloneCall = "gene+nt", 
                      method="ward.D2")
 
scatterClonotype(combined, 
                 cloneCall ="gene", 
                 x.axis = "t1_t", 
                 y.axis = "p1_p",
                 dot.size = "total",
                 graph = "proportion")

scatterClonotype(combined, 
                 cloneCall ="gene", 
                 x.axis = "t2_t", 
                 y.axis = "p2_p",
                 dot.size = "total",
                 graph = "proportion")

scatterClonotype(combined, 
                 cloneCall ="gene", 
                 x.axis = "t3_t", 
                 y.axis = "p3_p",
                 dot.size = "total",
                 graph = "proportion")

setwd("~/liuzhuqing")

a1 <- Read10X(data.dir = "~/liuzhuqing/Paratumor_1")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p1_p.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Paratumor_2")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p2_p.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Paratumor_3")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p3_p.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Tumor_1")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "t1_t.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Tumor_2")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "t2_t.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Tumor_3")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "t3_t.rds")

p1_p<-readRDS(file="p1_p.rds")
p2_p<-readRDS(file="p2_p.rds")
p3_p<-readRDS(file="p3_p.rds")
t1_t<-readRDS(file="t1_t.rds")
t2_t<-readRDS(file="t2_t.rds")
t3_t<-readRDS(file="t3_t.rds")

p1_p<-RenameCells(p1_p,add.cell.id="p1_p",for.merge=T)
p1_p@meta.data$tech<-"para"
p1_p@meta.data$celltype<-"para_p1"

p2_p<-RenameCells(p2_p,add.cell.id="p2_p",for.merge=T)
p2_p@meta.data$tech<-"para"
p2_p@meta.data$celltype<-"para_p2"

p3_p<-RenameCells(p3_p,add.cell.id="p3_p",for.merge=T)
p3_p@meta.data$tech<-"para"
p3_p@meta.data$celltype<-"para_p3"

t1_t<-RenameCells(t1_t,add.cell.id="t1_t",for.merge=T)
t1_t@meta.data$tech<-"tumor"
t1_t@meta.data$celltype<-"tumor_t1"

t2_t<-RenameCells(t2_t,add.cell.id="t2_t",for.merge=T)
t2_t@meta.data$tech<-"tumor"
t2_t@meta.data$celltype<-"tumor_t2"

t3_t<-RenameCells(t3_t,add.cell.id="t3_t",for.merge=T)
t3_t@meta.data$tech<-"tumor"
t3_t@meta.data$celltype<-"tumor_t3"

t12<-merge(t1_t,t2_t)
t123<-merge(t12,t3_t)
p12<-merge(p1_p,p2_p)
p123<-merge(p12,p3_p)
tp<-merge(t123,p123)

saveRDS(tp, file="tp_before_integrate.rds")

hms<-readRDS(file="tp_before_integrate.rds")

#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:4)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)

#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}

reference.list <- pancreas.list[c( "tumor_t1","tumor_t2", "tumor_t3",
                                   "para_p1","para_p2","para_p3")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:4)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")


hms_individual_integrated<-readRDS(file="hms_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1

#find how many 11 clusters
ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:10)
hms_cluster <- FindClusters( hms_neighbor, resolution = 0.5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:4)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test.rds")

hms_cluster<-readRDS(file="hms_cluster_test.rds")
cluster0.markers <- FindMarkers(hms_cluster, ident.1=0, min.pcr=0.25)
head(cluster0.markers, n=10)
cluster1.markers <- FindMarkers(hms_cluster, ident.1=1, min.pcr=0.25)
head(cluster1.markers, n=10)
cluster2.markers <- FindMarkers(hms_cluster, ident.1=2, min.pcr=0.25)
head(cluster2.markers, n=10)
cluster3.markers <- FindMarkers(hms_cluster, ident.1=3, min.pcr=0.25)
head(cluster3.markers, n=10)
cluster4.markers <- FindMarkers(hms_cluster, ident.1=4, min.pcr=0.25)
head(cluster4.markers, n=10)
cluster5.markers <- FindMarkers(hms_cluster, ident.1=5, min.pcr=0.25)
head(cluster5.markers, n=10)
cluster6.markers <- FindMarkers(hms_cluster, ident.1=6, min.pcr=0.25)
head(cluster6.markers, n=10)
cluster7.markers <- FindMarkers(hms_cluster, ident.1=7, min.pcr=0.25)
head(cluster7.markers, n=10)
cluster8.markers <- FindMarkers(hms_cluster, ident.1=8, min.pcr=0.25)
head(cluster8.markers, n=10)
cluster9.markers <- FindMarkers(hms_cluster, ident.1=9, min.pcr=0.25)
head(cluster9.markers, n=10)
cluster10.markers <- FindMarkers(hms_cluster, ident.1=10, min.pcr=0.25)
head(cluster10.markers, n=10)
cluster11.markers <- FindMarkers(hms_cluster, ident.1=11, min.pcr=0.25)
head(cluster11.markers, n=10)
cluster12.markers <- FindMarkers(hms_cluster, ident.1=12, min.pcr=0.25)
head(cluster12.markers, n=10)
cluster13.markers <- FindMarkers(hms_cluster, ident.1=13, min.pcr=0.25)
head(cluster13.markers, n=10)
cluster14.markers <- FindMarkers(hms_cluster, ident.1=14, min.pcr=0.25)
head(cluster14.markers, n=10)
cluster15.markers <- FindMarkers(hms_cluster, ident.1=15, min.pcr=0.25)
head(cluster15.markers, n=10)
cluster16.markers <- FindMarkers(hms_cluster, ident.1=16, min.pcr=0.25)
head(cluster16.markers, n=10)
cluster17.markers <- FindMarkers(hms_cluster, ident.1=17, min.pcr=0.25)
head(cluster17.markers, n=10)
cluster18.markers <- FindMarkers(hms_cluster, ident.1=18, min.pcr=0.25)
head(cluster18.markers, n=10)

new.cluster.ids <- c("Mesenchymal","T helper","B","Mesenchymal","Mesenchymal",
                      "Endothelial","CD8 T","Mesenchymal","Th17","Macrophage","Epithelial","Dendritic","Mesenchymal",
                     "NKT", "Endothelial","Basal","Monocyte","Treg","Treg") 
setwd("~/liuzhuqing")
hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")

#each type of cells
Mesenchymal<-subset(hms_cluster_id, idents=c('Mesenchymal'))
DimPlot(Mesenchymal, reduction = "umap")
saveRDS(Mesenchymal, file="Mesenchymal.rds")

T_helper<-subset(hms_cluster_id, idents=c('T helper'))
DimPlot(T_helper, reduction = "umap")
saveRDS(T_helper, file="T_helper.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Endothelial<-subset(hms_cluster_id, idents=c('Endothelial'))
DimPlot(Endothelial, reduction = "umap")
saveRDS(Endothelial, file="Endothelial.rds")

CD8_T<-subset(hms_cluster_id, idents=c('CD8 T'))
DimPlot(CD8_T, reduction = "umap")
saveRDS(CD8_T, file="CD8_T.rds")

Th17<-subset(hms_cluster_id, idents=c('Th17'))
DimPlot(Th17, reduction = "umap")
saveRDS(Th17, file="Th17.rds")

Macrophage<-subset(hms_cluster_id, idents=c('Macrophage'))
DimPlot(Macrophage, reduction = "umap")
saveRDS(Macrophage, file="Macrophage.rds")

Epithelail<-subset(hms_cluster_id, idents=c('Epithelial'))
DimPlot(Epithelail, reduction = "umap")
saveRDS(Epithelail, file="Epithelail.rds")

Dendritic<-subset(hms_cluster_id, idents=c('Dendritic'))
DimPlot(Dendritic, reduction = "umap")
saveRDS(Dendritic, file="Dendritic.rds")

NKT<-subset(hms_cluster_id, idents=c('NKT'))
DimPlot(NKT, reduction = "umap")
saveRDS(NKT, file="NKT.rds")

Basal<-subset(hms_cluster_id, idents=c('Basal'))
DimPlot(Basal, reduction = "umap")
saveRDS(Basal, file="Basal.rds")

Monocyte<-subset(hms_cluster_id, idents=c('Monocyte'))
DimPlot(Monocyte, reduction = "umap")
saveRDS(Monocyte, file="Monocyte.rds")

Treg<-subset(hms_cluster_id, idents=c('Treg'))
DimPlot(Treg, reduction = "umap")
saveRDS(Treg, file="Treg.rds")


names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

#add seurat
setwd("~/liuzhuqing")
seurat <- readRDS(file="hms_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))

setwd("~/liuzhuqing/tcr")

t1<- read.csv("t1_filtered_contig_annotations.csv")
t2<- read.csv("t2_filtered_contig_annotations.csv")
t3<- read.csv("t3_filtered_contig_annotations.csv")
p1<- read.csv("p1_filtered_contig_annotations.csv")
p2<- read.csv("p2_filtered_contig_annotations.csv")
p3<- read.csv("p3_filtered_contig_annotations.csv")

contig_list <- list(t1,t2,t3,p1,p2,p3)
combined <- combineTCR(contig_list, 
                       samples = c( "t1", "t2","t3","p1", "p2", "p3"), 
                       ID = c( "t","t","t","p", "p", "p" ), cells ="T-AB")

seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene", group.by = "sample", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

names(seurat@meta.data)
head(seurat@meta.data)

DimPlot(seurat, group.by = "tech") +
  scale_color_manual(values=colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())

DimPlot(seurat, group.by = "celltype") +
  scale_color_manual(values=colorblind_vector(8), na.value="grey") + 
  theme(plot.title = element_blank())

slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (100 < X <= 500)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))
DimPlot(seurat, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())

clonalOverlay(seurat, reduction = "umap", 
              freq.cutpoint = 30, bins = 10, facet = "tech") + 
  guides(color = "none")

clonalOverlay(seurat, reduction = "umap", 
              freq.cutpoint = 30, bins = 10, facet = "celltype") + 
  guides(color = "none")

#clonalNetwork
#ggraph needs to be loaded due to issues with ggplot
library(ggraph)

#No Identity filter
clonalNetwork(seurat, 
              reduction = "umap", 
              identity = "ident",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")

#Examining NKT only
clonalNetwork(seurat, 
              reduction = "umap", 
              identity = "ident",
              filter.identity = "NKT",
              cloneCall = "aa")

#Examining Treg only
clonalNetwork(seurat, 
              reduction = "umap", 
              identity = "ident",
              filter.identity = "Treg",
              cloneCall = "aa")

#Examining CD8_T only
clonalNetwork(seurat, 
              reduction = "umap", 
              identity = "ident",
              filter.identity = "CD8_T",
              cloneCall = "aa")

#highlightClonotypes

seurat <- highlightClonotypes(seurat, 
                              cloneCall= "aa", 
                              sequence = c("CAVRDGDYKLSF_CASSQDLAGGTDTQYF", 
                                           "CAGSDQGAQKLVF_CASSVGGGRYTGELFF",
                                           "CAASLFQGAQKLVF_CASSSDGVEAFF",
                                           "CAVRDGDYKLSF_CASSQDLAGGTDTQYF",
                                           "CACTDQTGANNLFF_CASSLGGTMNTEAFF",
                                           "CAESINQAGTALIF_CATSRASLRAQYF",
                                           "CAGGGSGAGSYQLTF_CASSFGLPNYGYTF",
                                           "CALSEARGETSGSRLTF_CASSMGQRTDTQYF",
                                           "CALSEARNSNYQLIW_CASSFDLTGSGNTIYF",
                                           "CALSRSGAGSYQLTF_CAGDREETQYF",
                                           "CALVGYSSASKIIF;CAVRDGDYKLSF_CASSQDLAGGTDTQYF",
                                           "CALVTQGGSEKLVF_CASSLSGTGELFF",
                                           "CAMNSGGYQKVTF_CATEDSNYGYTF",
                                           "CAMNSGGYQKVTF_NA",
                                           "CAMNSGGYQKVTF;CAVILPRGNNDMRF_CATEDSNYGYTF",
                                           "CAMRERDGVGLTGGGNKLTF_CATSRDELQREIAKNIQYF",
                                           "CATDAEKLVF_CASSLMGAGAADTQYF",
                                           "CATDILIEFSGGYNKLIF_CASSLSSGVGCYEQYF",
                                           "CAVEDPSTGANNLFF_CASSETGSSTDTQYF",
                                           "CAVEDYGGSQGNLIF_CASSLASGGNIQYF",
                                           "CAVNARLMF_CASSLPLGLAMNTEAFF",
                                           "CAYRIKAAGNKLTF_CASSGGSSYEQYF",
                                           "CIALYSGYALNF;CAVNARLMF_CASSLPLGLAMNTEAFF",
                                           "CIVRLAGTALIF_CASSSSPGLALDGELFF",
                                           "CIVRVEDMRF;CAYRIKAAGNKLTF_CASSGGSSYEQYF",
                                           "NA_CASSETGSSTDTQYF",
                                           "NA_CASSVGQPGLAGGETQYF"))
DimPlot(seurat, group.by = "highlight",pt.size = 0.4)
  theme(plot.title = element_blank())
  
  
NKT <- combineExpression(combined, NKT, 
                              cloneCall="gene", group.by = "sample", proportion = FALSE, 
                              cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))  

NKT<- highlightClonotypes(NKT, 
                                cloneCall= "aa", 
                                sequence = c("CAVRDGDYKLSF_CASSQDLAGGTDTQYF", 
                                             "CAGSDQGAQKLVF_CASSVGGGRYTGELFF",
                                             "CAASLFQGAQKLVF_CASSSDGVEAFF",
                                             "CAVRDGDYKLSF_CASSQDLAGGTDTQYF",
                                             "CACTDQTGANNLFF_CASSLGGTMNTEAFF",
                                             "CAESINQAGTALIF_CATSRASLRAQYF",
                                             "CAGGGSGAGSYQLTF_CASSFGLPNYGYTF",
                                             "CALSEARGETSGSRLTF_CASSMGQRTDTQYF",
                                             "CALSEARNSNYQLIW_CASSFDLTGSGNTIYF",
                                             "CALSRSGAGSYQLTF_CAGDREETQYF",
                                             "CALVGYSSASKIIF;CAVRDGDYKLSF_CASSQDLAGGTDTQYF",
                                             "CALVTQGGSEKLVF_CASSLSGTGELFF",
                                             "CAMNSGGYQKVTF_CATEDSNYGYTF",
                                             "CAMNSGGYQKVTF_NA",
                                             "CAMNSGGYQKVTF;CAVILPRGNNDMRF_CATEDSNYGYTF",
                                             "CAMRERDGVGLTGGGNKLTF_CATSRDELQREIAKNIQYF",
                                             "CATDAEKLVF_CASSLMGAGAADTQYF",
                                             "CATDILIEFSGGYNKLIF_CASSLSSGVGCYEQYF",
                                             "CAVEDPSTGANNLFF_CASSETGSSTDTQYF",
                                             "CAVEDYGGSQGNLIF_CASSLASGGNIQYF",
                                             "CAVNARLMF_CASSLPLGLAMNTEAFF",
                                             "CAYRIKAAGNKLTF_CASSGGSSYEQYF",
                                             "CIALYSGYALNF;CAVNARLMF_CASSLPLGLAMNTEAFF",
                                             "CIVRLAGTALIF_CASSSSPGLALDGELFF",
                                             "CIVRVEDMRF;CAYRIKAAGNKLTF_CASSGGSSYEQYF",
                                             "NA_CASSETGSSTDTQYF",
                                             "NA_CASSVGQPGLAGGETQYF"))
DimPlot(NKT01, group.by = "highlight",pt.size = 0.4)
  theme(plot.title = element_blank())
  
Treg <- combineExpression(combined, Treg, 
                           cloneCall="gene", group.by = "sample", proportion = FALSE, 
                           cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))  
  
Treg<- highlightClonotypes(Treg, 
                               cloneCall= "aa", 
                               sequence = c("CAVRDGDYKLSF_CASSQDLAGGTDTQYF", 
                                            "CAGSDQGAQKLVF_CASSVGGGRYTGELFF",
                                            "CAASLFQGAQKLVF_CASSSDGVEAFF",
                                            "CAVRDGDYKLSF_CASSQDLAGGTDTQYF",
                                            "CACTDQTGANNLFF_CASSLGGTMNTEAFF",
                                            "CAESINQAGTALIF_CATSRASLRAQYF",
                                            "CAGGGSGAGSYQLTF_CASSFGLPNYGYTF",
                                            "CALSEARGETSGSRLTF_CASSMGQRTDTQYF",
                                            "CALSEARNSNYQLIW_CASSFDLTGSGNTIYF",
                                            "CALSRSGAGSYQLTF_CAGDREETQYF",
                                            "CALVGYSSASKIIF;CAVRDGDYKLSF_CASSQDLAGGTDTQYF",
                                            "CALVTQGGSEKLVF_CASSLSGTGELFF",
                                            "CAMNSGGYQKVTF_CATEDSNYGYTF",
                                            "CAMNSGGYQKVTF_NA",
                                            "CAMNSGGYQKVTF;CAVILPRGNNDMRF_CATEDSNYGYTF",
                                            "CAMRERDGVGLTGGGNKLTF_CATSRDELQREIAKNIQYF",
                                            "CATDAEKLVF_CASSLMGAGAADTQYF",
                                            "CATDILIEFSGGYNKLIF_CASSLSSGVGCYEQYF",
                                            "CAVEDPSTGANNLFF_CASSETGSSTDTQYF",
                                            "CAVEDYGGSQGNLIF_CASSLASGGNIQYF",
                                            "CAVNARLMF_CASSLPLGLAMNTEAFF",
                                            "CAYRIKAAGNKLTF_CASSGGSSYEQYF",
                                            "CIALYSGYALNF;CAVNARLMF_CASSLPLGLAMNTEAFF",
                                            "CIVRLAGTALIF_CASSSSPGLALDGELFF",
                                            "CIVRVEDMRF;CAYRIKAAGNKLTF_CASSGGSSYEQYF",
                                            "NA_CASSETGSSTDTQYF",
                                            "NA_CASSVGQPGLAGGETQYF"))
DimPlot(Treg, group.by = "highlight",pt.size = 0.4)
theme(plot.title = element_blank())
  
CD8 <- combineExpression(combined, CD8_T, 
                          cloneCall="gene", group.by = "sample", proportion = FALSE, 
                          cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))  

CD8<- highlightClonotypes(CD8, 
                           cloneCall= "aa", 
                           sequence = c("CAVRDGDYKLSF_CASSQDLAGGTDTQYF", 
                                        "CAGSDQGAQKLVF_CASSVGGGRYTGELFF",
                                        "CAASLFQGAQKLVF_CASSSDGVEAFF",
                                        "CAVRDGDYKLSF_CASSQDLAGGTDTQYF",
                                        "CACTDQTGANNLFF_CASSLGGTMNTEAFF",
                                        "CAESINQAGTALIF_CATSRASLRAQYF",
                                        "CAGGGSGAGSYQLTF_CASSFGLPNYGYTF",
                                        "CALSEARGETSGSRLTF_CASSMGQRTDTQYF",
                                        "CALSEARNSNYQLIW_CASSFDLTGSGNTIYF",
                                        "CALSRSGAGSYQLTF_CAGDREETQYF",
                                        "CALVGYSSASKIIF;CAVRDGDYKLSF_CASSQDLAGGTDTQYF",
                                        "CALVTQGGSEKLVF_CASSLSGTGELFF",
                                        "CAMNSGGYQKVTF_CATEDSNYGYTF",
                                        "CAMNSGGYQKVTF_NA",
                                        "CAMNSGGYQKVTF;CAVILPRGNNDMRF_CATEDSNYGYTF",
                                        "CAMRERDGVGLTGGGNKLTF_CATSRDELQREIAKNIQYF",
                                        "CATDAEKLVF_CASSLMGAGAADTQYF",
                                        "CATDILIEFSGGYNKLIF_CASSLSSGVGCYEQYF",
                                        "CAVEDPSTGANNLFF_CASSETGSSTDTQYF",
                                        "CAVEDYGGSQGNLIF_CASSLASGGNIQYF",
                                        "CAVNARLMF_CASSLPLGLAMNTEAFF",
                                        "CAYRIKAAGNKLTF_CASSGGSSYEQYF",
                                        "CIALYSGYALNF;CAVNARLMF_CASSLPLGLAMNTEAFF",
                                        "CIVRLAGTALIF_CASSSSPGLALDGELFF",
                                        "CIVRVEDMRF;CAYRIKAAGNKLTF_CASSGGSSYEQYF",
                                        "NA_CASSETGSSTDTQYF",
                                        "NA_CASSVGQPGLAGGETQYF"))
DimPlot(CD8, group.by = "highlight",pt.size = 0.4)
theme(plot.title = element_blank())

  
#occupiedscRepertoire
occupiedscRepertoire(seurat, x.axis = "ident")
  
#alluvialClonotypes
alluvialClonotypes(seurat, cloneCall = "gene", 
                   y.axes = c("celltype", "ident", "tech"), 
                   color = "ident") 
  
#getCirclize
library(circlize)
library(scales)

circles <- getCirclize(seurat, 
                       group.by = "ident")

#Just assigning the normal colors to each cluster
grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)

#Graphing the chord diagram
circlize::chordDiagram(circles,
                       self.link = 1, 
                       grid.col = grid.cols)


subset <- subset(seurat, tech == "tumor")

circles <- getCirclize(subset, group.by = "ident")

grid.cols <- scales::hue_pal()(length(unique(subset@active.ident)))
names(grid.cols) <- levels(subset@active.ident)

chordDiagram(circles, self.link = 1, 
             grid.col = grid.cols, directional = 1, 
             direction.type =  "arrows",
             link.arr.type = "big.arrow")

subset <- subset(seurat, tech == "para")

circles <- getCirclize(subset, group.by = "ident")

grid.cols <- scales::hue_pal()(length(unique(subset@active.ident)))
names(grid.cols) <- levels(subset@active.ident)

chordDiagram(circles, self.link = 1, 
             grid.col = grid.cols, directional = 1, 
             direction.type =  "arrows",
             link.arr.type = "big.arrow")

#Clustering Clonotypes
sub_combined <- clusterTCR(combined[[2]], 
                           chain = "TRA", 
                           sequence = "aa", 
                           threshold = 0.85, 
                           group.by = NULL)

#TRB_cluster
seurat <- clusterTCR(seurat,
                     chain = "TRB",
                     group.by = "celltype", 
                     sequence = "aa", 
                     threshold = 0.85)

DimPlot(seurat, group.by = "TRB_cluster") + 
  scale_color_manual(values = colorblind_vector(length(unique(seurat@meta.data[,"TRB_cluster"])))) + 
  NoLegend()

#TRA_cluster
seurat <- clusterTCR(seurat,
                     chain = "TRA",
                     group.by = "celltype", 
                     sequence = "aa", 
                     threshold = 0.85)

DimPlot(seurat, group.by = "TRA_cluster") + 
  scale_color_manual(values = colorblind_vector(length(unique(seurat@meta.data[,"TRA_cluster"])))) + 
  NoLegend()


#StartracDiversity
StartracDiversity(seurat, 
                  type = "tech", 
                  sample = "celltype", 
                  by = "overall")

#Clonotype Bias  not work
clonotypeBias(seurat, 
              cloneCall = "aa", 
              split.by = "celltype", 
              group.by = "ident",
              n.boots = 20, 
              min.expand =10)

clonotypeBias(seurat, 
              cloneCall = "aa", 
              split.by = "Patient", 
              group.by = "ident",
              n.boots = 20, 
              min.expand =10)


#Error in dimnames(x) <- dn : 
#  length of 'dimnames' [2] not equal to array extent

#Working with clonotypes after clustering Clonal Diversity
combined2 <- expression2List(seurat, 
                             split.by = "ident")
length(combined2) #now listed by cluster


#Clonal Diversity
clonalDiversity(combined2, 
                cloneCall = "nt")

#Clonal Homeostasis
clonalHomeostasis(combined2, 
                  cloneCall = "nt")

#Clonal Proportion
clonalProportion(combined2, 
                 cloneCall = "nt")

#Clonal Overlap
clonalOverlap(combined2, 
              cloneCall="aa", 
              method="overlap")
#
clonalDiversity(seurat, 
                split.by = "celltype", 
                cloneCall = "nt")

#
clonalDiversity(seurat, 
                split.by = "tech", 
                cloneCall = "nt")

