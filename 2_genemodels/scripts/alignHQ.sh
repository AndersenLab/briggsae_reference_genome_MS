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

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"

#cd $wkdir/${file%%.*}
#mkdir alignHQ
cd $wkdir/${file%%.*}/alignments/PacBio

source activate cupcake
minimap2 -ax splice:hq -t 12 -uf \
   $wkdir/$file $wkdir/longreads/${file%%.*}/cluster/${file%%.*}.pol.hq.fasta > ${file%%.*}.hq.transcripts.sam

source activate bedtools
samtools sort -o ${file%%.*}.hq.transcripts.sorted.bam -O bam -@ 8 ${file%%.*}.hq.transcripts.sam
