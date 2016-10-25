#!/bin/bash

set -ex

usage() {
    echo "usage: submitRNASeqPreprocessing.sh replicate [T1, T2 or T3]
" 1>&2
    exit 1
}

if [ $# != 1 ]; then 
	usage
fi

runproj=b2010064

mail="david.sundell@umu.se"
in=/proj/b2015062/sequence_data/raw/$1
out=/proj/b2015062/nobackup/SpruceCrossSection
genome=/proj/b2011227/indices/STAR/Pabies01-genome
#gff3=/proj/b2011227/reference/gff3/Eugene.gff3
gff3=~/final_filter_other.gff3
start=8
end=8
mem=fat

export SORTMERNADIR=/home/davidsu/opt/sortmerna-2.0

module load bioinfo-tools FastQC samtools star htseq

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fastq.gz"`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
#newgrp b2013155 
bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -m $mem -g $genome -G $gff3 -H $gff3 -a $runproj $mail

'''
required:
		<proj>			project Name
		<mail>			user mail
		<dir>			Path to the PASA folder 
		<pasaDB>		Name of the PASA database
		<genome>		Path to the species genome fasta
	optional:
		-C 				Specify a database version for PASA default (latest)
		-i 				Add to look for only intergenic novel genes.
		-d 				dry run
		-c 				specify frameDP CFG file, default cfg in P.tricocarpa ref
		-t 				Train dataset for frameDP, in case no previous training has 
						been done on the species. takes path to the species fasta file
'''