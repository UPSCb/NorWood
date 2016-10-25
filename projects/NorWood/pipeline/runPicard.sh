#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 4
#SBATCH -t 0-12:00:00
#SBATCH --mail-type=ALL

# -p node is needed to accept the -C memory configuration

## stop on error and be verbose in the output
set -e -x

## load the modules
module load bioinfo-tools Picard-tools/2.4.1

for bam in /mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/STAR/all_bam/*STAR.bam
do

	basename=`basename $bam`
	java -jar /mnt/picea/Modules/apps/bioinfo/picard/2.4.1/picard.jar CollectInsertSizeMetrics I=$bam O=$1/insert_size_metrics$basename.txt H=$1/insert_size_histogram$basename.pdf
	exit
done
