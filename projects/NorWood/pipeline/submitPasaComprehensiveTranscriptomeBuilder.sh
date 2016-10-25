#!/bin/bash -l

## error and verbose
set -ex

## vars
mail="rhedensjo@gmail.com"

## define a function
usage () {
echo >&2 \
"This function uses a .cfg file
"
exit 1
}

## args number
if [ $# != 1 ]; then
    usage
    exit 1
fi

config=/mnt/picea/projects/u2015027/cfg/Pabies01.cfg
out=/mnt/picea/projects/u2015027/PASA/rna/`date +%Y%m%d`

## submit
sbatch --mail-user=$mail -e $out/pasa-build.err -o $out/pasa-build.out -J P-Build-$1 $UPSCb/pipeline/runPasaComprehensiveTranscriptomeBuilder.sh $config $out
