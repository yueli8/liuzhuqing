library(Seurat)

#volcano_plot
setwd("~/liuzhuqing/volcano_data")

res <- read.csv("cd8t_volcano.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-6,6),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)


res <- read.csv("nkt_volcano.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-6,6),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

res <- read.csv("Th17_volcano.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-8,7),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

res <- read.csv("Thelper_volcano.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-6,6),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

res <- read.csv("Treg_volcano.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-8,6),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

#deg of CD8,NKT,Th17,T_helper and Treg
setwd("~/liuzhuqing")
hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")
DimPlot(hms_cluster_id, label = T) + NoLegend()
DimPlot(hms_cluster_id, label = T) 
table(Idents(hms_cluster_id))

Mesenchymal<-subset(hms_cluster_id, idents=c('Mesenchymal'))
DimPlot(Mesenchymal, reduction = "umap")
saveRDS(Mesenchymal, file="Mesenchymal.rds")

Thelper<-subset(hms_cluster_id, idents=c('T_helper'))
DimPlot(Thelper, reduction = "umap")
saveRDS(Thelper, file="Thelper.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Endothelial<-subset(hms_cluster_id, idents=c('Endothelial'))
DimPlot(Endothelial, reduction = "umap")
saveRDS(Endothelial, file="Endothelial.rds")

CD8<-subset(hms_cluster_id, idents=c('CD8_T'))
DimPlot(CD8, reduction = "umap")
saveRDS(CD8, file="CD8.rds")

Th17<-subset(hms_cluster_id, idents=c('Th17'))
DimPlot(Th17, reduction = "umap")
saveRDS(Th17, file="Th17.rds")

Macrophage<-subset(hms_cluster_id, idents=c('Macrophage'))
DimPlot(Macrophage, reduction = "umap")
saveRDS(Macrophage, file="Macrophage.rds")

Epithelail<-subset(hms_cluster_id, idents=c('Epithelial'))
DimPlot(Epithelail, reduction = "umap")
saveRDS(Epithelail, file="Epithelial.rds")

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

#input each cluster
hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")
DimPlot(hms_cluster_id, label = T) + NoLegend()
DimPlot(hms_cluster_id, label = T) 
table(Idents(hms_cluster_id))

Mesenchymal<-readRDS("Mesenchymal.rds")
Thelper<-readRDS("Thelper.rds")
B<-readRDS("B.rds")
Endothelial<-readRDS("Endothelial.rds")
CD8<-readRDS("CD8.rds")
Th17<-readRDS("Th17.rds")
Macrophage<-readRDS("Macrophage.rds")
Epithelail<-readRDS("Epithelial.rds")
Dendritic<-readRDS("Dendritic.rds")
NKT<-readRDS("NKT.rds")
Basal<-readRDS("Basal.rds")
Monocyte<-readRDS("Monocyte.rds")
Treg<-readRDS("Treg.rds")

#deg in NKT
NKT<-readRDS("NKT.rds")
a<-NKT@meta.data
write.table(a,"a")
class(NKT)
NKT.sec<-as.SingleCellExperiment(NKT)
group<-factor(c(rep(1,177),rep(2,281)))

rds<-readRDS('NKT.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_NKT")
write.table(results.classified,"results.classified_NKT")

#deg in Treg
Treg<-readRDS("Treg.rds")
b<-Treg@meta.data
write.table(b,"b")
class(Treg)
Treg.sec<-as.SingleCellExperiment(Treg)
group<-factor(c(rep(1,128),rep(2,344)))

counts<-as.matrix(Treg@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_Treg")
write.table(results.classified,"results.classified_Treg")

#deg in CD8_T
CD8_T<-readRDS("CD8_T.rds")
c<-CD8_T@meta.data
write.table(c,"c")
class(CD8_T)
CD8_T.sec<-as.SingleCellExperiment(CD8_T)
group<-factor(c(rep(1,515),rep(2,465)))

counts<-as.matrix(CD8_T@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD8T")
write.table(results.classified,"results.classified_CD8T")

#deg in Thelper
Thelper<-readRDS("T_helper.rds")
e<-Thelper@meta.data
write.table(e,"e")
class(Thelper)
Thelper.sec<-as.SingleCellExperiment(Thelper)
group<-factor(c(rep(1,694),rep(2,908)))

counts<-as.matrix(Thelper@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_Thelper")
write.table(results.classified,"results.classified_Thelper")

#deg in Th17
Th17<-readRDS("Th17.rds")
d<-Th17@meta.data
write.table(d,"d")
class(Th17)
Th17.sec<-as.SingleCellExperiment(Th17)
group<-factor(c(rep(1,582),rep(2,226)))

counts<-as.matrix(Th17@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_Th17")
write.table(results.classified,"results.classified_Th17")


