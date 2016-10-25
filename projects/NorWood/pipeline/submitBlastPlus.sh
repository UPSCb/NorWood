#!/bin/bash -l

## error and verbose
set -e
set -x

## check for the UPSCb env. var.                                                                                                  
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
    usage
fi

## default args
#in=/mnt/picea/projects/spruce/Sitka_FLcDNA/fasta
#out=/mnt/picea/projects/spruce/Sitka_FLcDNA/blast/PacBio
#inxDir=/mnt/picea/storage/reference/Picea-abies/v1.1/GeneSpace/indices/blast
#inx=Pabies-PacBio-cDNA.fa

### split this file Potri03_update.12239_features.fasta end with fna
in=/mnt/picea/projects/aspseq/wood/PASA/20150211/fasta/
#out=/mnt/picea/projects/aspseq/wood/blast/Potra01
#inxDir=/mnt/picea/storage/reference/Populus-tremula/v1.0/indices/BLAST+/
#inx=Potra01-genome-with-header.fa

## Blast against NCBI
out=/mnt/picea/projects/aspseq/wood/blast/NCBI
inxDir=/mnt/picea/storage/reference/NCBI/20150318
#inx=nt

blastcommand="blastp"

if [ ! $1 == "" ]; then
	in=$1
fi
if [ ! $2 == "" ]; then
	inxDir=$2
fi
if [ ! $3 == "" ]; then
	out=$3
fi
if [ ! $4 == "" ]; then
	inx=$4
fi
if [ ! $5 == "" ]; then
	blastcommand=$5
fi


## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
for file in `find $in -maxdepth 1 -name "*.fasta"`; do
	sbatch -n 4 -e $out/`basename ${file//.fasta/.err}` -o $out/`basename ${file//.fasta/.out}` $UPSCb/pipeline/runBlastPlus.sh -p 4 $blastcommand $file $inxDir/$inx $out
done

