#!/bin/bash

## set
set -ex

## global vars
mail='rhedensjo@gmail.com'

## job var
gff=/mnt/picea/storage/reference/Picea-abies/v1.0/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene.gff3

## submit job
sbatch --mem=2G --mail-user $mail -e ${gff3//.gff/-validation.err} -o ${gff3//.gff/-validation.out} -J validate-`basename $gff` $UPSCb/pipeline/runPASA_GFF3_validator.sh $gff