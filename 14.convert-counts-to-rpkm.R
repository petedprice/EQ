library("edgeR")

data <- read.table("read_counts.txt",stringsAsFactors=F,header=T, row.names=1)
#data FORMAT: Geneid	Sample1	Sample2
#data FORMAT: MSTRG.43175	145.0	285.0
names(data)
dim(data)
conditions <- factor(c())
#conditions eg conditions <- factor(c("M","M","M","M","M","M","M","M","M","M","F","F","F","F","M"))
length(conditions)

#Extract RPKM
expr <- DGEList(counts=data)
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
rpkm <- rpkm(expr, log=FALSE,gene.length=gene_length_vector)
str(rpkm)
write.table(rpkm, file="rpkm",quote=F, sep="\t")

