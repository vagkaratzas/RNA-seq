#!/usr/bin/env Rscript
#setwd("C:/Users/vagos/Desktop/Bio/R/rnaseq")

library(QuasR)
sampleFile1 <- "samples_fastq.txt"
sampleFile2 <- "samples_bam.txt"
genomeFile <- "BSgenome.Hsapiens.NCBI.GRCh38"
proj1 <- qAlign(sampleFile1, genomeFile) #qAlign(sampleFile1, genomeFile, splicedAlignement=TRUE) #using SpliceMap instead of bowtie
write.table(alignments(proj1)$genome, sampleFile2, sep="/t", row.names=FALSE)
proj2 <- qAlign(sampleFile2, genomeFile, paired='no')
#####
qQCReport(proj2, "qc_report.pdf") #fastQC equivalent

library(rtracklayer)
library(GenomicFeatures)
annotFile <- "Homo_sapiens.GRCh38.83.chr.gtf"
#export(BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38, "BSgenome.Hsapiens.NCBI.GRCh38.fasta") #do once
chrLen <- scanFaIndex("BSgenome.Hsapiens.NCBI.GRCh38.fasta")
chrominfo <- data.frame(chrom=as.character(seqnames(chrLen)), length=width(chrLen), is_circular=rep(FALSE, length(chrLen)))
txdb <- makeTxDbFromGFF(file=annotFile, format="gtf", chrominfo=chrominfo, dataSource="NCBI", organism="Homo sapiens")
geneLevels <- qCount(proj2, txdb, reportLevel="gene") #check how txdb is made
#geneLevels #DEBUG
write.table(geneLevels, "gene_count.txt", row.names=TRUE, col.names=TRUE, sep="/t")
