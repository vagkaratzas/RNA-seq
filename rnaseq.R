#!/usr/bin/env Rscript

library(QuasR)
sampleFile1 <- "samples_fastq.txt"
sampleFile2 <- "samples_bam.txt"
genomeFile <- "BSgenome.Hsapiens.NCBI.GRCh38.fasta" #or BSgenome.Hsapiens.NCBI.GRCh38 to automaticaly download the genome
proj1 <- qAlign(sampleFile1, genomeFile) #qAlign(sampleFile1, genomeFile, splicedAlignement=TRUE) #uses SpliceMap instead of bowtie
write.table(alignments(proj1)$genome, sampleFile2, sep="\t", row.names=FALSE)
proj2 <- qAlign(sampleFile2, genomeFile, paired='no') #alternative qProject creation method, from bam files
#####
qQCReport(proj2, "qc_report.pdf") #fastQC equivalent

library(rtracklayer)
library(GenomicFeatures)
library(Rsamtools)
annotFile <- "Homo_sapiens.GRCh38.83.chr.gtf" #important: must be the same version as the alignment genome
#export(BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38, "BSgenome.Hsapiens.NCBI.GRCh38.fasta") #do once if you don't already have the .fasta genome
chrLen <- scanFaIndex("BSgenome.Hsapiens.NCBI.GRCh38.fasta")
chrominfo <- data.frame(chrom=as.character(seqnames(chrLen)), length=width(chrLen), is_circular=rep(FALSE, length(chrLen)))
txdb <- makeTxDbFromGFF(file=annotFile, format="gtf", chrominfo=chrominfo, dataSource="NCBI", organism="Homo sapiens")
geneLevels <- qCount(proj2, txdb, reportLevel="gene") #check how txdb is made
#geneLevels #DEBUG
write.table(geneLevels, "gene_count.txt", row.names=TRUE, col.names=TRUE, sep="\t")
