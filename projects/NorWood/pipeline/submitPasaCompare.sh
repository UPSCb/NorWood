#!/bin/bash -l

## error and verbose
set -ex

## vars
mail="rhedensjo@gmail.com"

## define a function
usage () {
echo >&2 \
"This function takes a .cfg file and a .gff3 file

"
exit 1
}

## args number
#if [ $# != 1 ]; then
#    usage
#    exit 1
#fi

## in case of a desired restart
restart=

## define vars
config=/mnt/picea/projects/u2015027/cfg/Pabies01.cfg
gff3=/mnt/picea/storage/reference/Picea-abies/v1.0/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene.gff3
out=/mnt/picea/projects/u2015027/PASA/rna/

## check
if [ ! -z $restart ];then
    restart="-r $restart"
fi
#Was before UPSCb/pipeline and after pasa-comp.out
#-J P-Comp-$1

## submit
sbatch --mail-user=$mail -e $out/pasa-comp.err -o $out/pasa-comp.out $UPSCb/pipeline/runPasaCompare.sh $restart $config $gff3 $out
