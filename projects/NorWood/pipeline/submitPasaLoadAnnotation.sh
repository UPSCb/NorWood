#!/bin/bash -l

## error and verbose
set -ex

## vars
mail="rhedensjo@gmail.com"

## define a function
usage () {
echo >&2 \
"This function uses .fa, .gff3, .cfg files
"
exit 1
}

## args number
if [ $# != 1 ]; then
    usage
    exit 1
fi

## global vars
#Ptremula
#fasta=/mnt/picea/storage/reference/Populus-tremula/v1.0/fasta/Potra01-genome.fa
#gff3=/mnt/picea/storage/reference/Populus-tremula/v1.0/gff3/Potra01-genome-forPasa-withDuplicatedExons.gff3

#Pabies-abies
fasta=/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenomeAssemblies/Pabies01-genome.fa
gff3=/mnt/picea/storage/reference/Picea-abies/v1.0/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene.gff3
config=/mnt/picea/projects/u2015027/cfg/Pabies01.cfg
out=/mnt/picea/projects/u2015027/PASA/rna/`date +%Y%m%d`


## submit
sbatch --mail-user=$mail -e $out/pasa-load.err -o $out/pasa-load.out -J P-Load-$1 $UPSCb/pipeline/runPasaLoadAnnotation.sh $config $fasta $gff3
