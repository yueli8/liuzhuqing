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

#export table, no graph
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
clonalDiversity(combined, 
                cloneCall = "gene", 
                group.by = "sample", 
                x.axis = "ID", 
                n.boots = 100)
 
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