#!/bin/bash -l

## be verbose and print
set -ex

## define a function
usage () {
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
}

## process the argument
in=/mnt/picea/projects/u2015027/trimmomatic
out=/mnt/picea/projects/u2015027/Trinity/`date +%Y%m%d`

## check vars
if [ -z $UPSCb ]; then
    echo "The UPSCb var needs to be set."
    usage
fi

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## find does not return _1 and _2 in the same order...
fwd=`find $in -name '*_trimmomatic_1.fq.gz' -printf %p\, | sed "s/,$//g"`
rev=`echo $fwd | sed "s/_1.fq.gz/_2.fq.gz/g"`

## execute
sbatch -w watson --mem=380G --mail-user="rhedensjo@gmail.com" -e $out/trinity.err -o $out/trinity.out -p core -n 48 -t 7-00:00:00 $UPSCb/pipeline/runTrinity.sh -n -p 48 -m 350G $out $fwd $rev

