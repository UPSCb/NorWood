#!/bin/bash

## global vars
run_proj=b2010064
proj=b2015062
mail="david.sundell@umu.se"
in=/proj/$proj/nobackup/SpruceCrossSection/$1
out=/proj/$proj/nobackup/SpruceCrossSection/$1
echo $in
## raw
for f in `find $in -name "*.fq.gz" -type f`; do
	fnam=`basename $f`
	sbatch -A $run_proj --mail-user $mail -e $out/fq_validator_${fnam//.fq.gz/.err} -o $out/fq_validator_${fnam//.fq.gz/.out} -J fqv-$fnam ../../../pipeline/runFastQValidator.sh $f
done





