#!/bin/bash
for i in *.bam;
do 
	cufflinks $i
	mv transcripts.gtf $i.transcripts.gtf
	mv isoforms.fpkm_tracking $i.isoforms.fpkm_tracking
	mv genes.fpkm_tracking $i.genes.fpkm_tracking
	mv skipped.gtf $i.skipped.gtf
done
