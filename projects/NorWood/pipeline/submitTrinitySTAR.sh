#!/bin/bash -l

## define a function
usage () {
echo >&2 \
"This function take two argument as parameters:
     the genome; use 'Pabies01'
     the trinity assembly; 'wood'
Note:You need to set the UPSCb env. variable to your UPSCb git checkout directory."
exit 1
}

## args number
if [ $# != 2 ]; then
    usage
fi

## check for the UPSCb env. var.
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
    usage
fi

## set the option
mail="rhedensjo@gmail.com"

case "$1" in
    Pabies01)
	genome=/mnt/picea/storage/reference/Picea-abies/v1.1/indices/STAR/Pabies01-genome/
	gtf=/mnt/picea/storage/reference/Picea-abies/v1.0/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene.gff3
	;;
    *) usage;;
esac


fwd=/mnt/picea/projects/u2015027/trinity/left.norm.fq_ext_all_reads.normalized_K25_C50_pctSD200.fq.gz
rev=/mnt/picea/projects/u2015027/trinity/right.norm.fq_ext_all_reads.normalized_K25_C50_pctSD200.fq.gz
out=/mnt/picea/projects/u2015027/STAR/Pabies01/
q
## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

sbatch -w watson --mem=256G --mail-user $mail -e $out/TrinitySTAR.err -o $out/TrinitySTAR.out -J STAR-T-$2 $UPSCb/pipeline/runSTAR.sh -q -g $gtf -f gff3 -l 30000000000 $out $genome $fwd $rev