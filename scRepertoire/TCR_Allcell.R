library(Seurat)
library(scRepertoire)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(dplyr)
library(Matrix)
#library(muscat)
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

#add seurat
setwd("~/liuzhuqing/herong")
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

#Examining T only
clonalNetwork(seurat, 
              reduction = "umap", 
              identity = "ident",
              filter.identity = "T",
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
#sub_combined <- clusterTCR(combined[[2]], 
#                       chain = "TRA", 
#                       sequence = "aa", 
#                       threshold = 0.85, 
#                       group.by = NULL)

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
new.cluster.ids <- c("Endothelial","Myofibroblast","Inflammatory_Fibroblast","Melanocyte",
                     "Inflammatory_Fibroblast","Endothelial","T","Inflammatory_Fibroblast",
                     "Myofibroblast","T","Dendritic","Myofibroblast","Inflammatory_Fibroblast",
                     "Melanocyte", "Inflammatory_Fibroblast","Epithelial","Epithelial",
                     "Natural_Killer","Dendritic","Melanocyte", "Inflammatory_Fibroblast",
                     "Endothelial", "Epithelial","Melanocyte", "Mast","Neutrophil",
                     "Epithelial","Endothelial", "Melanocyte","Inflammatory_Fibroblast",
                     "Inflammatory_Fibroblast", "Melanocyte","Epithelial")
names(new.cluster.ids) <- levels(seurat)
seurat <- RenameIdents(seurat, new.cluster.ids)
names(seurat@meta.data)
#[1] "orig.ident"             "nCount_RNA"             "nFeature_RNA"           "percent.mt"             "tech"                  
#[6] "celltype"               "integrated_snn_res.0.5" "seurat_clusters"        "barcode"                "CTgene"                
#[11] "CTnt"                   "CTaa"                   "CTstrict"               "Frequency"              "cloneType"   
head(seurat[[]])

clonotypeBias(seurat, 
              cloneCall = "aa", 
              split.by = "tech", 
              group.by = "RNA_snn_res.0.8",
              n.boots = 20, 
              min.expand =10)

clonotypeBias(seurat, 
              cloneCall = "aa", 
              split.by = "tech", 
              group.by = "seurat_clusters",
              n.boots = 20, 
              min.expand =10)

clonotypeBias(seurat, 
              cloneCall = "aa", 
              split.by = "tech", 
              group.by = "ident",
              n.boots = 20, 
              min.expand =10)


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

