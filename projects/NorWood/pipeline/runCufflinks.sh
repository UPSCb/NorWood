#!/bin/bash -l
#SBATCH -p core
#SBATCH -t 5-00:00:00
#SBATCH -n 16
#SBATCH --mem 40G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.sundell@umu.se

## stop on error and be verbose
set -ex

ncores=16
output_folder=$1
gtf_reference=$2
fa_file=$3
infile=$4

module load bioinfo-tools
module load cufflinks/2.2.1
#module load samtools/0.1.19

#mkdir "tmp"
#samtools view $infile | awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C$

cufflinks -o $output_folder -p $ncores -g $gtf_reference --library-type fr-unstranded -I 15000 --no-faux-reads -b $fa_file -F Potra_cuff $infile 
