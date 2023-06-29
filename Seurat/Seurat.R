library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(dplyr)
library(Matrix)
#library(muscat)
library(reshape2)
library(celldex)
library(clustree)
setwd("~/liuzhuqing/herong")

a1 <- Read10X(data.dir = "~/liuzhuqing/Paratumor_1")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "p1_p.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Paratumor_2")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "p2_p.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Paratumor_3")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "p3_p.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Tumor_1")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "t1_t.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Tumor_2")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "t2_t.rds")

a1 <- Read10X(data.dir = "~/liuzhuqing/Tumor_3")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:20)
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
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:20)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) 
p1+p2

#pbmc <- JackStraw(pancreas, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)

ElbowPlot(pancreas)
hms_neighbor<- FindNeighbors(pancreas, dims = 1:20)
#obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
#clustree(obj)
hms_cluster <- FindClusters(hms_neighbor, resolution = 0.8)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:20)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "my_hms_cluster_test.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.25, 
                                logfc.threshold = 0.25
)
write.table(scRNA.markers,file="mycellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="mytop20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="mytop20_marker_genes.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(hms)
#assay
hms@assays
dim(hms@meta.data)
View(hms@meta.data)

#我的結果是0-30簇,何蓉的結果是0-32簇,我應該在她的文件上標注.
setwd("~/liuzhuqing/herong")
hms_cluster<-readRDS(file="herong_hms_cluster_test.rds")
DimPlot(hms_cluster, reduction = "umap")

new.cluster.ids <- c("Endothelial","Myofibroblast","Inflammatory_Fibroblast","Melanocyte",
                     "Inflammatory_Fibroblast","Endothelial","T","Inflammatory_Fibroblast",
                     "Myofibroblast","T","Dendritic","Myofibroblast","Inflammatory_Fibroblast",
                     "Melanocyte", "Inflammatory_Fibroblast","Epithelial","Epithelial",
                     "Natural_Killer","Dendritic","Melanocyte", "Inflammatory_Fibroblast",
                     "Endothelial", "Epithelial","Melanocyte", "Mast","Neutrophil",
                     "Epithelial","Endothelial", "Melanocyte","Inflammatory_Fibroblast",
                     "Inflammatory_Fibroblast", "Melanocyte","Epithelial") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
names(hms_cluster_id@meta.data)

DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

setwd("~/liuzhuqing/herong")
hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))

#each type of cells
Endothelial<-subset(hms_cluster_id, idents=c('Endothelial'))
DimPlot(Endothelial, reduction = "umap")
saveRDS(Endothelial, file="Endothelial.rds")

myCAF<-subset(hms_cluster_id, idents=c('Myofibroblast'))
DimPlot(myCAF, reduction = "umap")
saveRDS(myCAF, file="myCAF.rds")

iCAF<-subset(hms_cluster_id, idents=c('Inflammatory_Fibroblast'))
DimPlot(iCAF, reduction = "umap")
saveRDS(iCAF, file="iCAF.rds")

Melanocyte<-subset(hms_cluster_id, idents=c('Melanocyte'))
DimPlot(Melanocyte, reduction = "umap")
saveRDS(Melanocyte, file="Melanocyte.rds")

T<-subset(hms_cluster_id, idents=c('T'))
DimPlot(T, reduction = "umap")
saveRDS(T, file="T.rds")

Dendritic<-subset(hms_cluster_id, idents=c('Dendritic'))
DimPlot(Dendritic, reduction = "umap")
saveRDS(Dendritic, file="Dendritic.rds")

Epithelial<-subset(hms_cluster_id, idents=c('Epithelial'))
DimPlot(Epithelial, reduction = "umap")
saveRDS(Epithelial, file="Epithelial.rds")

NK<-subset(hms_cluster_id, idents=c('Natural_Killer'))
DimPlot(NK, reduction = "umap")
saveRDS(NK, file="NK.rds")

Mast<-subset(hms_cluster_id, idents=c('Mast'))
DimPlot(Mast, reduction = "umap")
saveRDS(Mast, file="Mast.rds")

Neutrophil<-subset(hms_cluster_id, idents=c('Neutrophil'))
DimPlot(Neutrophil, reduction = "umap")
saveRDS(Neutrophil, file="Neutrophil.rds")

#readRDS
Endothelial<-readRDS(file="Endothelial.rds")
myCAF<-readRDS(file="myCAF.rds")
iCAF<-readRDS(file="iCAF.rds")
Melanocyte<-readRDS(file="Melanocyte.rds")
T<-readRDS(file="T.rds")
Dendritic<-readRDS(file="Dendritic.rds")
Epithelial<-readRDS(file="Epithelial.rds")
NK<-readRDS(file="NK.rds")
Mast<-readRDS(file="Mast.rds")
Neutrophil<-readRDS(file="Neutrophil.rds")

#regroup T cell
T<-readRDS(file="T.rds")
DimPlot(T, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dim(T@meta.data)
T <- NormalizeData(T, normalization.method = "LogNormalize", scale.factor = 10000)
T <- FindVariableFeatures(T, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(T)
#T- ScaleData(T, rownames(T))#跑不出來
T <- RunPCA(T, features = VariableFeatures(object = T))
#Determine the ‘dimensionality’ of the dataset
#T <- JackStraw(T, num.replicate = 100,dims = 40)
#T <- ScoreJackStraw(T, dims = 1:40)
#JackStrawPlot(T, dims = 1:40)
#ElbowPlot(T,ndims = 40)
#確定下面的dims
hms_neighbor<- FindNeighbors(T, dims = 1:10)
#obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#clustree(obj)
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
hms_cluster <- FindClusters(hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:10)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "myT_cluster_test.rds")

hms_cluster<-readRDS(file="myT_cluster_test.rds")
DimPlot(hms_cluster, reduction = "umap")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.1, #如果設定是0.25,第11個聚類只有14個基因,無法write.csv
                                logfc.threshold = 0.25
)
write.table(scRNA.markers,file="myT_cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="myT_top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(top20_table,"myT_top20_marker_genes.csv",row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(T)
#assay
T@assays
dim(T@meta.data)
View(T@meta.data)

#T cell annotation,use herong's data.
setwd("~/liuzhuqing/herong")
hms_cluster<-readRDS(file="herong_T_cluster_test .rds")
DimPlot(hms_cluster, reduction = "umap")

new.cluster.ids <- c("T_helper","CD4_T","Treg","CD8_T","Naive_T","CD8_T","MAIT",
                     "Naive_T","CD8_T","NKT","Treg","CD8_T","CD8_T",
                     "Melanocyte","Endothelial","B") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
names(hms_cluster_id@meta.data)

DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "herongT_cluster_id_test.rds")

setwd("~/liuzhuqing/herong")
hms_cluster_id<-readRDS(file="herongT_cluster_id_test.rds")

#each type of cells
T_helper<-subset(hms_cluster_id, idents=c('T_helper'))
DimPlot(T_helper, reduction = "umap")
saveRDS(T_helper, file="T_helper.rds")

CD4_T<-subset(hms_cluster_id, idents=c('CD4_T'))
DimPlot(CD4_T, reduction = "umap")
saveRDS(CD4_T, file="CD4_T.rds")

Treg<-subset(hms_cluster_id, idents=c('Treg'))
DimPlot(Treg, reduction = "umap")
saveRDS(Treg, file="Treg.rds")

CD8_T<-subset(hms_cluster_id, idents=c('CD8_T'))
DimPlot(CD8_T, reduction = "umap")
saveRDS(CD8_T, file="CD8_T.rds")

Naive_T<-subset(hms_cluster_id, idents=c('Naive_T'))
DimPlot(Naive_T, reduction = "umap")
saveRDS(Naive_T, file="Naive_T.rds")

MAIT<-subset(hms_cluster_id, idents=c('MAIT'))
DimPlot(MAIT, reduction = "umap")
saveRDS(MAIT, file="MAIT.rds")

NKT<-subset(hms_cluster_id, idents=c('NKT'))
DimPlot(NKT, reduction = "umap")
saveRDS(NKT, file="NKT.rds")

Melanocyte<-subset(hms_cluster_id, idents=c('Melanocyte'))
DimPlot(Melanocyte, reduction = "umap")
saveRDS(Melanocyte, file="Melanocyte.rds")

Endothelial<-subset(hms_cluster_id, idents=c('Endothelial'))
DimPlot(Endothelial, reduction = "umap")
saveRDS(Endothelial, file="Endothelial.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

#violin plot_Tmarker
VlnPlot(hms_cluster_id, features = c("CD8A", "CD8B","GZMK","CD27","CCR7","IL7R","CD4","FOXP3"))
VlnPlot(hms_cluster_id, features = c("CTLA4","TRDC","TGFB1","TRAV1-2","LAG3","PDCD1","TGFBR2",
                                     "IGHG1","MZB1"))
#ridgeplot
RidgePlot(T, features =c("TTN"))

#feature plot
FeaturePlot(T, features = c("GBP5"))

#新增代码：去除非T细胞
T_new <- subset(hms_cluster_id, idents = c("T_helper","CD4_T","Treg","CD8_T",
                                    "Naive_T","MAIT","NKT"))
saveRDS(T_new, file = "t.rds")
DimPlot(T_new, reduction = "umap", label = FALSE, pt.size = 0.5) 
DimPlot(T_new, reduction = "umap", label = TRUE, pt.size =0.5) 

#each type of cells
T_helper<-subset(T_new, idents=c('T_helper'))
DimPlot(T_helper, reduction = "umap")
saveRDS(T_helper, file="T_helper.rds")

CD4_T<-subset(T_new, idents=c('CD4_T'))
DimPlot(CD4_T, reduction = "umap")
saveRDS(CD4_T, file="CD4_T.rds")

Treg<-subset(T_new, idents=c('Treg'))
DimPlot(Treg, reduction = "umap")
saveRDS(Treg, file="Treg.rds")

CD8_T<-subset(T_new, idents=c('CD8_T'))
DimPlot(CD8_T, reduction = "umap")
saveRDS(CD8_T, file="CD8_T.rds")

Naive_T<-subset(T_new, idents=c('Naive_T'))
DimPlot(Naive_T, reduction = "umap")
saveRDS(Naive_T, file="Naive_T.rds")

MAIT<-subset(T_new, idents=c('MAIT'))
DimPlot(MAIT, reduction = "umap")
saveRDS(MAIT, file="MAIT.rds")

NKT<-subset(T_new, idents=c('NKT'))
DimPlot(NKT, reduction = "umap")
saveRDS(NKT, file="NKT.rds")

setwd("~/liuzhuqing/herong")
hms_cluster<-readRDS(file="t.rds")
DimPlot(hms_cluster, reduction = "umap")

VlnPlot(hms_cluster, features = c("CD8A", "CD8B","GZMK","CD27","CCR7","IL7R","CD4","FOXP3"))
VlnPlot(hms_cluster, features = c("CTLA4","TRDC","TGFB1","TRAV1-2","LAG3","PDCD1","TGFBR2",
                                     "IGHG1","MZB1"))

RidgePlot(hms_cluster, features = c("CD8A", "CD8B","GZMK","CD27","CCR7","IL7R","CD4","FOXP3"))
RidgePlot(hms_cluster, features = c("CTLA4","TRDC","TGFB1","TRAV1-2","LAG3","PDCD1","TGFBR2",
                                  "IGHG1","MZB1"))

#feature plot
FeaturePlot(hms_cluster, features = c("CD8A", "CD8B","GZMK","CD27","CCR7","IL7R","CD4","FOXP3"))
FeaturePlot(hms_cluster, features = c("CTLA4","TRDC","TGFB1","TRAV1-2","LAG3","PDCD1","TGFBR2",
                                    "IGHG1","MZB1"))
