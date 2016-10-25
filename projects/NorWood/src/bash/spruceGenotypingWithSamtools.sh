#!/bin/bash -l
#SBATCH -p node
#SBATCH -c 1
#SBATCH -t 0-12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.sundell@umu.se
#SBATCH -A u2015027

module load bioinfo-tools bcftools samtools htslib

genome="/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenomeAssemblies/Pabies01-genome.fa"
outdir="/mnt/picea/home/david/genotyping/spruce/mergegenotypes"
indir="/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/STAR/bam"
refs="/mnt/picea/projects/u2015027/genotyping/reference/1_140729_BC49TKACXX_P1169_105_sortmerna_trimmomatic_STAR.bam /mnt/picea/projects/u2015027/genotyping/reference/1_140729_BC49TKACXX_P1169_107_sortmerna_trimmomatic_STAR.bam /mnt/picea/projects/u2015027/genotyping/reference/1_140729_BC49TKACXX_P1169_106_sortmerna_trimmomatic_STAR.bam"

fin=$1
#samples="T1-06_sortmerna_trimmomatic_STAR.bam "$indir"T2-06_sortmerna_trimmomatic_STAR.bam T3-06_sortmerna_trimmomatic_STAR.bam"
FILES=$indir/T[1,2,3]-06*STAR.bam
samples=""
for f in $FILES
do
	#echo $f
	samples="$samples $f"
done
#echo $samples

mergedVCF=$outdir/merged_withref.vcf
if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi
if [ ! -f $mergedVCF ]; then
	cat "" > $mergedVCF
fi
echo "Output file: "$mergedVCF
while read scaffold; do
  samtools mpileup -d100000 -gf $genome -r "$scaffold"$samples $refs | bcftools call -v -c - > $outdir/$scaffold.vcf
  grep -v "#" $outdir/$scaffold.vcf >> $mergedVCF
done < $fin