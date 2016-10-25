#!/bin/bash -l

proj=b2015062
tmp=/glob/davidsu/tmp
mail="david.sundell@umu.se"
out=/proj/$proj/nobackup/test/sortmerna
in=/proj/$proj/sequence_data/raw/

if [ ! -d $out ]; then
    mkdir -p $out
fi

export SORTMERNADIR=/home/davidsu/opt/sortmerna

for f in `find $in -name "*_[1,2].fastq.gz" -type l`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do 
fnam=`basename $line`
sbatch -A $proj --mail-user $mail -C "usage_mail" -e $out/$fnam.err -o $out/$fnam.out -J smr-$fnam ../../../pipeline/runSortmerna.sh $out $tmp ${line}_1.fastq.gz ${line}_2.fastq.gz
exit
done

