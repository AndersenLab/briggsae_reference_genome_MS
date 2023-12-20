#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="stringtie"

source activate stringtie

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"

cd $wkdir/${file%%.*}/predictions/stringtie/

stringtie $wkdir/${file%%.*}/alignments/PacBio/${file%%.*}.hq.transcripts.sorted.bam -o ${file%%.*}.stringtie.gtf -v -p 12
