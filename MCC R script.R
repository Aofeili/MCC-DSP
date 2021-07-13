

# Run the first time to install:
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")

BiocManager::install(c("ComplexHeatmap"))

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
BiocManager::install(
        
#old--------------------------------------------------
install.packages("heatmap.plus")


filename <- "Decont2.txt"
getwd()

# Read the data into a data.frame
my_data <- read.table(filename, sep="\t", quote="",
                      row.names=1, 
                      stringsAsFactors=FALSE,header=TRUE)


# Default parameters
Heatmap((my_data), 
        #cluster_columns = FALSE, 
        row_names_gp=gpar(cex=0.2),
        row_names_side = "left",
        column_names_side= "top",
        column_names_gp=gpar(cex=0.2 #, col= c(rep("brown3",10), rep("blue",15), rep("purple",17), rep("cyan4", 8)) 
                             ),
        width = unit(6,"cm"),
        height = unit(80,"cm"),
        #clustering_distance_rows ="maximum",
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D"
        )
?heatmap

#-PCA---------------------------------------------------------------------------

if (!requireNamespace('BiocManager', quietly = TRUE))
        install.packages('BiocManager')

BiocManager::install('PCAtools')
library(PCAtools)

p <- pca(my_data, metadata = NULL, removeVar = 0.1
         )
?pca
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p)
biplot(p, showLoadings = TRUE,
       labSize = 2, pointSize = 5, sizeLoadingsNames = 2)

library(cluster.datasets)
library(pca3d)
data(my_data)
my_data

pca <- prcomp(my_data[,-1], scale.=TRUE)
gr <- factor(my_data[-1,])
summary(gr)
head(my_data)
pca3d(pca, group=gr)

pca2d(pca, group=gr, legend="topleft", )


#-example--------------------------------------------------------
BiocManager::install('preprocessCore')
library(preprocessCore)
getwd()
filename <- "Example.txt"
data_f <- read.table(filename, sep="\t", quote="",
                     row.names=1, 
                     stringsAsFactors=FALSE,header=TRUE)
data_f
data_f_lg <- log(data_f)
data_f_lg


#data--------------------------------------------------
library(preprocessCore)
getwd()

filename <- "CK skin.txt"
data_f <- read.table(filename, sep="\t", quote="",
                      row.names=1, 
                      stringsAsFactors=FALSE,header=TRUE)

filename <- "CK pheno noskin.txt"
pheno_f <- read.table(filename, sep="\t", quote="",
                     row.names=1, 
                     stringsAsFactors=FALSE,header=TRUE)
pheno_f
#data_f
#rownames(data_f)
boxplot(data_f)
?boxplot

data_m<-as.matrix(data_f)
pheno_m<-as.matrix(pheno_f)

post.norm <- normalize.quantiles(data_m)
rownames(post.norm)<-rownames(data_f)
colnames(post.norm)<-colnames(data_f)

post.norm_f<-as.data.frame(post.norm)
#post.norm_f
boxplot(post.norm_f)

post.norm_lg <- log(post.norm)
write.csv(post.norm_lg, file="quantile normed.csv")

library(Biobase)
#eset<-ExpressionSet(assayData=post.norm_f,
 #                   phenoData=AnnotatedDataFrame(pheno_f),
  #                  featrueData=AnnotatedDataFrame(pheno_f))
#design<- model.matrix(~inflame, data=subset(pheno_f,tumor=="mcc"))

#pheno_sub=subset(pheno_f,block==101287)
#data_sub<-post.norm_lg[rownames(pheno_sub),]
#data_sub
#rownames(pheno_sub)
head(post.norm_lg)

design<- model.matrix(~inflame, data=pheno_sub)
                      
design<- model.matrix(~inflame, data=pheno_f)
design<- model.matrix(~tumor, data=pheno_f)
                                            
design<- model.matrix(~progression, data=pheno_f)
design
head(design)
#olSums(design)
#able(pheno_f[,"inflame"])
BiocManager::install('limma')
library(limma)
fit<-lmFit(post.norm,design)
fit<-eBayes(fit)
fit
results<-decideTests(fit[,"progressionp"])
results<-decideTests(fit[,"inflamett"])
results<-decideTests(fit[,"tumorskin"])

summary(results)
results
topTable_export<- topTable(fit, 
                           number=300,
                           resort.by="logFC",
                           coef=ncol(design))
write.csv(topTable_export, file="export7.csv")
?topTable


head(topTable_export)
idx=rownames(topTable_export)
idx=c("MMP3","CCL19","DNAJC14","ESR2","ACVR1C","THBS4","RXRG","IL19","PGM2","TDO2","BAMBI","MCAT","ASPG","APOA1","PEBP1","SOX10","CTSG","MIF","TTC30A","PYCR2","HACD2","ANP32B","CCR3","MLF1","ERCC2","FOXM1","RAD21","HPRT1"
)
  
TopGenes<-post.norm_lg[idx,]
TopGenes
dim(TopGenes)


head(post.norm_lg)
head(TopGenes)

?write.csv
library(ComplexHeatmap)
Heatmap((post.norm_lg), #quantile normed data
        #cluster_columns = FALSE, 
        row_names_gp=gpar(cex=0.2),
        row_names_side = "left",
        column_names_side= "top",
        column_names_gp=gpar(cex=0.2 #, col= c(rep("brown3",10), rep("blue",15), rep("purple",17), rep("cyan4", 8)) 
        ),
        width = unit(6,"cm"),
        height = unit(80,"cm"),
        #clustering_distance_rows ="maximum",
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D"
)

Heatmap((TopGenes), #topgenes
        cluster_columns = FALSE, 
        cluster_rows =  FALSE, 
        
        #row_split = rep(c("A", "B"), 9)
        row_names_gp=gpar(cex=0.2),
        row_names_side = "left",
        column_names_side= "top",
        column_names_gp=gpar(cex=0.3 #, col= c(rep("brown3",10), rep("blue",15), rep("purple",17), rep("cyan4", 8)) 
        ),
        width = unit(6,"cm"),
        height = unit(25,"cm"),
        clustering_distance_rows ="maximum",
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D"
)
?heatmap
Heatmap((data_m), #pre-quantile normed example
        #cluster_columns = FALSE, 
        #row_names_gp=gpar(cex=0.2),
        row_names_side = "left",
        column_names_side= "top",
        #column_names_gp=gpar(cex=0.2 #, col= c(rep("brown3",10), rep("blue",15), rep("purple",17), rep("cyan4", 8)) 
        #),
        #width = unit(6,"cm"),
        #height = unit(80,"cm"),
        #clustering_distance_rows ="maximum",
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D"
)
