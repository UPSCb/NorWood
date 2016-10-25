#!/bin/bash -l

## vars
proj=b2015062
mail="david.sundell@umu.se"

## define a function
usage(){
    echo "This function take one argument as parameter; one of 'raw','trimmomatic','sortmerna'"
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
    exit 1
}

## args number
if [ $# != 1 ]; then
    	echo "Error: wrong number of arguments"
	usage
fi

## process the argument
dir=
case "$1" in
    raw)
	dir="/proj/$proj/nobackup/SpruceCrossSection/fastqc/raw"
	;;
    trimmomatic) 
	dir="/proj/$proj/nobackup/SpruceCrossSection/fastqc/trimmomatic"
	;;
    sortmerna)
	dir="/proj/$proj/nobackup/SpruceCrossSection/fastqc/sortmerna"
	;;
esac

## stop if no dir
if [ -z $dir ]; then
    echo "Error: no directory supplied"
    usage
fi

## check that the dir exists
if [ ! -d $dir ]; then
    echo "Error: the directory supplied doest exist"
    usage
fi

## check that the UPSCb env var exists
if [ -z $UPSCb ]; then
    echo "Error: the UPSCb env var doesnt exist"
    usage
fi

## submit
bash $UPSCb/pipeline/runFastQCMultiviewer.sh $dir
