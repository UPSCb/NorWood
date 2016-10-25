#!/bin/bash -l

## be verbose and print
set -ex

## define a function
usage () {
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
}

## process the argument
fw=/mnt/picea/projects/u2015027/trinity/left.norm.fq_ext_all_reads.normalized_K25_C50_pctSD200.fq.gz
rv=/mnt/picea/projects/u2015027/trinity/right.norm.fq_ext_all_reads.normalized_K25_C50_pctSD200.fq.gz
out=/mnt/picea/projects/u2015027/trinity/GenomeGuided/trinity_`date +%Y%m%d`
genome=/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenomeAssemblies/Pabies01-genome.fa
bam=/mnt/picea/projects/u2015027/STAR/Pabies01/left.norm.fq_ext_all_reads.normalized_K25_C50_pctSD200.fq.gz_STAR.bam

user="rhedensjo"

## check vars
if [ -z $UPSCb ]; then
    echo "The UPSCb var needs to be set."
    usage
fi

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## execute
sbatch -w watson --mem=400G -n 48 --mail-user=$user"@gmail.com" -e $out/trinityGG.err -o $out/trinityGG.out -t 7-00:00:00 $UPSCb/pipeline/runTrinityGenomeGuided.sh -p 48 -m 240 $out $fw $rv $genome $bam