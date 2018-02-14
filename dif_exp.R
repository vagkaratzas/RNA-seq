#!/usr/bin/env Rscript
#setwd("C:/Users/vagos/Desktop/Bio/R/rnaseq")
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")

#counts <- read.table("gene_count.txt", header = TRUE, row.names = 1, sep = "\t")
#counts_matrix <- matrix(unlist(counts), ncol = 5, byrow = TRUE)
#rownames(counts_matrix) <- rownames(counts)
#colnames(counts_matrix) <- colnames(counts)
#counts_matrix2 <- counts_matrix[,-1] #removing width
x <- read.delim("gene_count.txt",row.names="Symbol")
x <- x[,-1] #removing width
library(edgeR)
group <- factor(c(1,2,1,2)) #doesn't work with only 1 pattern for each group, needs repetitions
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design) #estimates dispersion
# To perform quasi-likelihood F-tests:
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
tT <- topTags(qlf, n=100, sort.by="logFC")
#To perform likelihood ratio tests:
#fit <- glmFit(y,design)
#lrt <- glmLRT(fit,coef=2)
#topTags(lrt)
write.table(tT, "gene_dif_exp.txt", row.names=TRUE, col.names=TRUE, sep="\t")
