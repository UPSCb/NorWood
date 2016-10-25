#!/bin/bash -l

## error and verbose
set -ex

## module
## module load myPipe

## define a function
usage () {
echo >&2 \
"This function take one mandatory argument as parameter:
     the index to use; one of 'Pabies01'
A second argument (the gmap output format) can be optionally provided, e.g. samse. It defaults to gff3_gene.
Note:You need to set the UPSCb env. variable to your UPSCb git checkout directory."
exit 1
}

## args number
if [ "$#" -lt "1" ] || [ "$#" -gt "2" ]; then
    usage
    exit 1
fi

fmt="gff3_gene"
if [ $# -eq 2 ]; then
    fmt=$2
fi

## check for the UPSCb env. var.
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
    usage
fi

## default args
in=/mnt/picea/projects/u2015027/trinity/Trinity.fasta.gz
out=/mnt/picea/projects/u2015027/gmap

case "$1" in
Pabies1.0)
inxDir=/mnt/picea/storage/reference/Populus-abies/v1.0/indices/GMAP/$1
;;
esac


## create the out dir
if [ ! -d $out/$1 ]; then
    mkdir -p $out/$1
fi

## prepare
sbatch -p core -e $out/$1/$1.err -o $out/$1/$1.out $UPSCb/pipeline/runGMAP.sh -f $fmt $in $inxDir $1 $out/$1