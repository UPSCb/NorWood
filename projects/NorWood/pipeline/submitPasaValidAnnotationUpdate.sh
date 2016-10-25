#!/bin/bash -l

## error and verbose
set -ex

## vars
mail="rhedensjo@gmail.com"

## define a function
usage () {
echo >&2 \
"This function ueses a genome file 
"
exit 1
}

## args number
if [ $# != 1 ]; then
    usage
    exit 1
fi

options="-V -R"
## Include an output for failed sequeneces
options="-A"
options=""

genome=/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenomeAssemblies/Pabies01-genome.fa
out=/mnt/picea/projects/aspseq/asp201/rna/PASA/`date +%Y%m%d`
#	options="-s 10 $options"


## prepare
if [ ! -d $out ]; then
    mkdir -p $out
fi

echo sbatch --mail-user=$mail -e $out/pasa-dump_v2.err -o $out/pasa-dump_v2.out -J Pasa-dump-$1 $UPSCb/pipeline/runPasaValidAnnotationUpdate.sh $dbname $out $options

# /mnt/picea/home/david/Git/UPSCb/pipeline/runPasaValidAnnotationUpdate.sh Potri03 /mnt/picea/projects/aspseq/wood/PASA/20150211
