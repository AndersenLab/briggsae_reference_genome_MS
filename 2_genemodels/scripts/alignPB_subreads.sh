#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="alignHQ"

wkdir="/projects/b1059/projects/Nicolas/Clay"

#cd $wkdir/${file%%.*}
#mkdir alignHQ
#cd alignHQ

source activate cupcake
cd $wkdir/${file%%.*}/alignments/DNA/
minimap2 -ax map-pb $wkdir/genomes/$file $wkdir/${file%%.*}/alignments/DNA/${file%%.*}_all_reads.fastq > ${file%%.*}.filtered.subreads.sam

#source activate bedtools
#samtools sort -o ${file%%.*}.hq.transcripts.sorted.bam -O bam -@ 8 ${file%%.*}.hq.transcripts.sam
