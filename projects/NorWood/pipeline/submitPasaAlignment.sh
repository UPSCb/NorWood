#!/bin/bash -l

## error and verbose
set -ex

## vars
mail="rhedensjo@gmail.com"

## define a function
usage () {
echo >&2 \
"This function runs PASA Alingnment and uses .fa, .cfg, .gtf, Trinity.fasta, Trinity-GG.fasta files
"
exit 1
}

## global vars
fasta=/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenomeAssemblies/Pabies01-genome-collapsed-for-STAR.fa
options="-c 48 -s 23"


config=/mnt/picea/projects/u2015027/cfg/Pabies01.cfg
cufflinks=/mnt/picea/projects/u2015027/cufflink/left.norm.fq_ext_all_reads.normalized_K25_C50_pctSD200.fq.gz_STAR/transcripts.gtf
trinity=/mnt/picea/projects/u2015027/trinity/Trinity.fasta.gz
trinityGG=/mnt/picea/projects/u2015027/trinity/GenomeGuided/trinity_20160623/Trinity-GG.fasta.gz
out=/mnt/picea/projects/u2015027/PASA/rna/20160629

## prepare
if [ ! -d $out ]; then
    mkdir -p $out
fi

sbatch --mail-user=$mail -e $out/pasa.err -n 48 -o $out/pasa.out -J P-Aln-$1 $UPSCb/pipeline/runPasaAlignment.sh $options $config $fasta $trinity $trinityGG $cufflinks $out
