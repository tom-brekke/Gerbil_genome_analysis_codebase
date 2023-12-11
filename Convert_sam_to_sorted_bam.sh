#!/bin/bash

#this takes a sam file, converts it to a bam file and sorts the bam file, finally it deletes the intermediates (a sam and a non-sorted bam)
#needs a sam file at the command line ($1)
#it can take a second input "q" which is the low bound mapping quality but this isn't necessary
				
sam=$1 # here is the sam input
if [ ! -z $2 ]
then	
	q=$2
	bam=${sam/.sam/.filt$2.bam} #this is the non-sorted bam
	sorted_bam=${bam/filt$2.bam/filt$2.sorted.bam} #this is the sorted bam
else
	q="0"
	bam=${sam/.sam/.tmp.bam} #this is the non-sorted bam
	sorted_bam=${bam/.tmp.bam/.sorted.bam} #this is the sorted bam
fi

#sorted_prefix=${sorted_bam/.bam/} #and for samtools sort, it just wants the prefix

run_cmd.sh "samtools view -b -S -q $q $sam > $bam" $bam
run_cmd.sh "samtools sort $bam -o $sorted_bam" $sorted_bam
run_cmd.sh "samtools index $sorted_bam" $sorted_bam.bai 
rm $sam #this removes the original sam file
rm $bam #and the non-sorted bam file.

#these were all samtools version 1.1, but that's different on Musculus so try this...