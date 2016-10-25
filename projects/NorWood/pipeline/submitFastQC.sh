#!/bin/bash

## vars
proj=b2015062
mail="david.sundell@umu.se"

## define a function
usage () {
    echo "This function take two arguments as parameter; one of 'raw','trimmomatic','sortmerna' and T1 T2 or T3"
}

## args number
if [ $# != 2 ]; then
    usage
    exit 1
fi

## process the argument
dir=
pattern=
typ=
case "$1" in
    raw)
	dir="/proj/$proj/sequence_data/$1/$2"
	pattern="*.fastq.gz"
	typ="l"
	;;
    trimmomatic) 
	dir="/proj/$proj/nobackup/normTest/$2"
	pattern="*_trimmomatic_[1,2].fq.gz"
	typ="f"
	;;
    sortmerna)
	dir="/proj/$proj/nobackup/normTest/$2"
	pattern="*_sortmerna_[1,2].fq.gz"
	typ="f"
	;;
esac

## stop if no dir
if [ -z $dir ]; then
    usage
    exit 1
fi

## check that the dir exists
if [ ! -d $dir ]; then
    usage
    exit 1
fi

## create the args
out=/proj/$proj/SpruceCrossSection/analysis/FastQC/$2

if [ ! -d $out ]; then
    mkdir -p $out
fi

for file in `find $dir -name $pattern -type $typ`; do 
fnam=`basename $file`
fnam=${fnam//.f*q.gz/}
sbatch -A $proj --mail-user $mail -C "usage_mail" -e $out/$fnam.err -o $out/$fnam.out -J FastQC-$fnam ../../../pipeline/runFastQC.sh $out $file
done
