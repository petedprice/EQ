library("edgeR")

data <- read.table("read_counts.txt",stringsAsFactors=F,header=T, row.names=1)
#data FORMAT: Geneid	Sample1	Sample2
#data FORMAT: MSTRG.43175	145.0	285.0
names(data)
dim(data)
conditions <- factor(c())
#conditions eg conditions <- factor(c("M","M","M","M","M","M","M","M","M","M","F","F","F","F","M"))
length(conditions)

#Check raw read count data
expr <- DGEList(counts=data,group=conditions)
plotMDS(expr)
plotSmear(expr,pair=NULL, de.tags=NULL, xlab="Average logCPM", ylab="logFC", pch=19, cex=0.2, smearWidth=0.5, panel.first=grid(), smooth.scatter=FALSE, lowess=FALSE)
cpm_expr <- cpm(expr)
sample1 <- density(log2(cpm_expr[,1]))
sample2 <- density(log2(cpm_expr[,2]))
plot(sample1, xlab = "CPM", ylab = "Density",type="l",lwd=2,main="Raw log2 cpm")
lines(sample2, type="l",lwd=2)

#Check normalised read count data
expr <- DGEList(counts=data,group=conditions)
norm_expr <- calcNormFactors(expr)
plotMDS(norm_expr)
plotSmear(norm_expr,pair=NULL, de.tags=NULL, xlab="Average logCPM", ylab="logFC", pch=19, cex=0.2, smearWidth=0.5, panel.first=grid(), smooth.scatter=FALSE, lowess=FALSE)
cpm_norm_expr <- cpm(norm_expr)
sample1 <- density(log2(cpm_norm_expr[,1]))
sample2 <- density(log2(cpm_norm_expr[,2]))
plot(sample1, xlab = "CPM", ylab = "Density",type="l",lwd=2,main="Raw log2 cpm")
lines(sample2, type="l",lwd=2)

#Normalise and extract rpkm
expr <- DGEList(counts=data)
norm_expr <- calcNormFactors(expr)
gene_length <- read.table("read_counts.genelength",stringsAsFactors=F)
names(gene_length)
dim(gene_length)
expressed_genes <- rownames(data)
length(expressed_genes)
gene_length <- subset(gene_length, V1 %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$V1),]
gene_length_vector <- c(gene_length$V2)
all(gene_length$V1 == rownames(expr))
#should print TRUE
rpkm_norm <- rpkm(norm_expr, log=FALSE,gene.length=gene_length_vector)
write.table(rpkm_norm, file="rpkm.normalised",quote=F, sep="\t")

#Heat maps
library(pheatmap)
library(pvclust)
palette2 <-colorRamps::"matlab.like2"(n=200)
bootstraps = pvclust(log2(rpkm_norm+1), method.hclust="average", method.dist="euclidean")
plot(bootstraps)
pheatmap(log2(rpkm_norm+1), show_colnames=T, show_rownames=F, color = palette2,clustering_distance_cols = "euclidean", clustering_method="average") 
