#!/bin/bash

#This script runs abundance_estimates_to_matrix from trinoate
#You neeed to specify the condition and tissue
# e.g sh abundance.sh GUT HM OV
#

cond_1=$1
cond_2=$2
cond_3=$3


#echo $cond_1$tissue"1_transcripts"

abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map none $cond_1.1.spiro_rsem.dir.genes.results $cond_1.2.spiro_rsem.dir.genes.results $cond_1.3.spiro_rsem.dir.genes.results $cond_2.1.spiro_rsem.dir.genes.results $cond_2.2.spiro_rsem.dir.genes.results $cond_2.3.spiro_rsem.dir.genes.results $cond_3.1.spiro_rsem.dir.genes.results $cond_3.2.spiro_rsem.dir.genes.results $cond_3.3.spiro_rsem.dir.genes.results 
