#install.packages(c("devtools","pkgload"))
#devtools::install_github("immunomind/immunarch",ref="dev")#sometime not work, because of internet.
#work on gnv at home computer,mm:1.
library(immunarch)
setwd("~/liuzhuqing/immunarch")
file_path="~/liuzhuqing/immunarch/data"
immdata<-repLoad(file_path)

#ensembl_list<-read.table("tmp01",header=FALSE)
a<-immdata$meta#meta data
write.csv(a,"a")
b<-top(immdata$data[[1]])
write.csv(b,"b")

#each file
p1_para<-immdata$data[[1]]
write.csv(p1_para,"p1_para")
p2_para<-immdata$data[[2]]
write.csv(p2_para,"p2_para")
p3_para<-immdata$data[[3]]
write.csv(p3_para,"p3_para")
t1_tumor<-immdata$data[[4]]
write.csv(t1_tumor,"t1_tumor")
t2_tumor<-immdata$data[[5]]
write.csv(t2_tumor,"t2_tumor")
t3_tumor<-immdata$data[[6]]
write.csv(t3_tumor,"t3_tumor")

#Calculate and visualise basic statistics

#length  # Visualise the length distribution of CDR3
repExplore(immdata$data, "lens") %>% vis()  
exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
p1 <- vis(exp_len)
p1
p4 <- vis(exp_len, .by = "Status", .meta = immdata$meta)
p4

# sample diversity 
# Visualise the Chao1 diversity of repertoires, grouped by the patient status
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta) 
repDiversity(immdata$data) %>% vis(.by = c("Status", "Patient"), .meta = immdata$meta)  

# number of clonotypes
exp_vol <- repExplore(immdata$data, .method = "volume")
p1 <- vis(exp_vol, .by = c("Status"), .meta = immdata$meta)
p2 <- vis(exp_vol, .by = c("Status", "Patient"), .meta = immdata$meta)
p1 + p2

exp_vol <- repExplore(immdata$data, .method = "volume")
by_vec <- c("P", "P", "P", "T", "T", "T")
p <- vis(exp_vol, .by = by_vec)
p

exp_cnt <- repExplore(immdata$data, .method = "count") 
p2 <- vis(exp_cnt)
p2

#Clonality analysis of immune repertoires
imm_pr <- repClonality(immdata$data, .method = "clonal.prop")
imm_pr
imm_pr %>% vis()

imm_hom <- repClonality(immdata$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
imm_hom
vis(imm_hom) + vis(imm_hom, .by = c("Status", "Patient"), .meta = immdata$meta)

imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top
vis(imm_top) + vis(imm_top, .by = "Status", .meta = immdata$meta)

imm_rare <- repClonality(immdata$data, .method = "rare")
imm_rare
vis(imm_rare) + vis(imm_rare, .by = "Status", .meta = immdata$meta)

#overlap
imm_ov1 <- repOverlap(immdata$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(immdata$data, .method = "morisita", .verbose = F)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)
p1 + p2
vis(imm_ov1, "circos")

p1 <- vis(imm_ov2, .text.size = 2.5, .signif.digits = 1)
p2 <- vis(imm_ov2, .text.size = 2, .signif.digits = 2)
p1 + p2

vis(imm_ov1, "heatmap2")
# Apply different analysis algorithms to the matrix of public clonotypes:
# "mds" - Multi-dimensional Scaling
#repOverlapAnalysis(imm_ov1, "mds")
#repOverlapAnalysis(imm_ov1, "mds") %>% vis()

repOverlapAnalysis(imm_ov1, "tsne") %>% vis()

#基因使用分析
gene_stats()
geneUsage(immdata$data[[1]]) %>% vis()  # Visualise the V-gene distribution for the first repertoire
imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T)
head(imm_gu)
vis(imm_gu, .title = "Gene usage", .text.size = 2)
imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = F)
head(imm_gu_js)
imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F)
head(imm_gu_cor)
p1 <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 3)
p2 <- vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 3)
p1
p2

vis(geneUsageAnalysis ( imm_gu , "cosine" , .verbose  =  F ))
vis(geneUsageAnalysis(imm_gu, "js+dbscan", .verbose = F))
imm_cl_pca <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .verbose = F)
imm_cl_mds <- geneUsageAnalysis(imm_gu, "js+mds+kmeans", .verbose = F)
imm_cl_tsne <- geneUsageAnalysis(imm_gu, "js+tsne+kmeans", .perp = .01, .verbose = F)
## Perplexity should be lower than K!
p1 <- vis(imm_cl_pca, .plot = "clust")
p2 <- vis(imm_cl_mds, .plot = "clust")
p3 <- vis(imm_cl_tsne, .plot = "clust")
p1
p2
p3

imm_cl_pca2 <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .k = 3, .verbose = F)
vis(imm_cl_pca2)

p1 <- vis(spectratype(immdata$data[[1]], .quant = "id", .col = "nt"))
p2 <- vis(spectratype(immdata$data[[1]], .quant = "count", .col = "aa+v"))
p1+p2
p2

# Compute statistics and visualise them
# Chao1 diversity measure
div_chao <- repDiversity(immdata$data, "chao1")
head(div_chao)
div_hill <- repDiversity(immdata$data, "hill")
head(div_hill)
div_d50 <- repDiversity(immdata$data, "d50")
head(div_d50)
div_div <- repDiversity(immdata$data, "div")
head(div_div)
p1 <- vis(div_chao)
p2 <- vis(div_chao, .by = c("Status", "Patient"), .meta = immdata$meta)
p3 <- vis(div_hill, .by = c("Status", "Patient"), .meta = immdata$meta)

p4 <- vis(div_d50)
p5 <- vis(div_d50, .by = "Status", .meta = immdata$meta)
p6 <- vis(div_div)
p1 + p2
p3 + p6
p4 + p5
 
repDiversity(immdata$data, "raref", .verbose = F) %>% vis(.log = TRUE)

tc1 <- trackClonotypes(immdata$data, list(1, 10), .col = "aa+v")
head(tc1)
p1  <-  vis ( tc1 )
p1

#only blood samples
file_path="~/liuzhuqing/immunarch/tumor"
immdata<-repLoad(file_path)
tc2  <-  trackClonotypes ( immdata$data , list (1,10), .col  =  "nt" )
head(tc2)
p2  <-  vis ( tc2 )
p2
tc2  <-  trackClonotypes ( immdata$data , list (1,10), .col  =  "aa+v" )
head(tc2)
p2  <-  vis ( tc2 )
p2

file_path="~/liuzhuqing/immunarch/paracarcinoma"
immdata<-repLoad(file_path)
tc2  <-  trackClonotypes ( immdata$data , list (1,10), .col  =  "nt" )
head(tc2)
p2  <-  vis ( tc2 )
p2
tc2  <-  trackClonotypes ( immdata$data , list (1,10), .col  =  "aa+v" )
head(tc2)
p2  <-  vis ( tc2 )
p2

file_path="~/liuzhuqing/immunarch/data"
immdata<-repLoad(file_path)
target <- immdata$data[[1]] %>%
select(CDR3.aa, V.name) %>%
head(10)
tc <- trackClonotypes(immdata$data, target)
vis(tc)

tc <- trackClonotypes(immdata$data, target, .col = "aa")
vis(tc, .plot = "smooth")
vis(tc, .plot = "area")
vis(tc, .plot = "line")

# Passing indices
names(immdata$data)[c(1, 3, 5)] # check sample names
vis(tc, .order = c(1, 3, 5))

immdata$meta$Timepoint <- sample(1:length(immdata$data))
immdata$meta

vdjdb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb")
vdjdb
#我们可以通过设置.species，.chain以及.pathology等参数进行过滤筛选出需要的信息
vdjdb  =  dbLoad ( "https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz" , "vdjdb" , .species  =  "HomoSapiens" , .chain  =  "TRB" , .pathology  =  "CMV" )
vdjdb
vdjdb_st = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/SearchTable-2019-10-17%2012_36_11.989.tsv.gz", "vdjdb-search", .species = "HomoSapiens", .chain = "TRB", .pathology = "CMV")
vdjdb_st
mcpas = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/McPAS-TCR.csv.gz", "mcpas", .species = "Human", .chain = "TRB", .pathology = "Cytomegalovirus (CMV)")
mcpas
#tbadb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/TBAdb.xlsx", "tbadb", .species = "Homo Sapiens", .chain = c("TRB", "TRA-TRB"), .pathology = "CMV")
#tbadb

dbAnnotate(immdata$data, vdjdb, "CDR3.aa", "cdr3")
dbAnnotate ( immdata $ data , mcpas , c ( "CDR3.aa" , "V.name" ), c ( "CDR3.beta.aa" , "TRBV" ))

#main command
repClonality(immdata$data, "homeo") %>% vis() 
repOverlap(immdata$data) %>% vis()  
geneUsage(immdata$data[[1]]) %>% vis()
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta) 

#kmers and motif
kmers  <-  getKmers ( immdata $ data [[ 1 ]], 5 )
kp  <-  kmer_profile ( kmers , "self" )
p1  <-  vis ( kp )
p2  <-  vis ( kp , .plot  =  "seq" )
p1 + p2
