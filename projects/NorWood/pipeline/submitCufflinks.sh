#!/bin/bash

## vars
proj=u2015027
out=/mnt/picea/projects/u2015027/cufflink/`date +%Y%m%d`
#inp=/mnt/picea/projects/u2015027/STAR/Pabies01/left.norm.fq_ext_all_reads.normalized_K25_C50_pctSD200.fq.gz_STAR.bam
inp=/mnt/picea/projects/u2015027/STAR/Pabies01/AlignedOut_noS.bam
fa=/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenomeAssemblies/Pabies01-genome-collapsed-for-STAR.fa
gtf=/mnt/picea/storage/reference/Picea-abies/v1.0/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene.gff3
mail="rhedensjo@gmail.com"
## arg
job="Pabies-abies"
## check
if [ -z $job ]; then
    echo "The argument should be Pabies-abies"
    exit 1
fi
if [ ! -d "$out" ]; then
	mkdir $out
fi
#for f in `find $inp -name "$pattern" -type f`
#do 
#	echo $f
#	fnam=`basename ${f//.bam/}`
#	 sbatch -A $proj --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out -J $job-$fnam ../../../pipeline/runCufflinks.sh -i 70000 $out $fa $inp $gtf
#done
sbatch -A $proj --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out -J $job-$fnam ../../../pipeline/runCufflinks.sh -i 70000 $out $fa $inp $gtf