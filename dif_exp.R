#!/usr/bin/env Rscript

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")

x <- read.delim("gene_count.txt",row.names="Symbol") #must have given the row name SYmbol manually inside gene counts
x <- x[,-1] #removing width
library(edgeR)
group <- factor(c(1,1,2,2)) #doesn't work with only 1 pattern for each group, needs repetitions
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design) #estimates dispersion
# To perform quasi-likelihood F-tests:
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
tT <- topTags(qlf, n=3000, sort.by="logFC")
#To perform likelihood ratio tests:
#fit <- glmFit(y,design)
#lrt <- glmLRT(fit,coef=2)
#topTags(lrt)
write.table(tT, "gene_dif_exp.txt", row.names=TRUE, col.names=TRUE, sep="\t")
