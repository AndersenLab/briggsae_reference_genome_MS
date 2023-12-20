#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="repmask"

source activate rep_mask

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"

cd $wkdir/${GENOME%%.*}/repeatmasker/soft
cp $wkdir/$GENOME $wkdir/${GENOME%%.*}/repeatmasker/soft/
RepeatMasker -s -xsmall -lib $wkdir/${GENOME%%.*}/blastx/${GENOME%%.*}.clust.class.noprot.fa -gff -pa 16 $wkdir/${GENOME%%.*}/repeatmasker/soft/$GENOME

cd $wkdir/${GENOME%%.*}/repeatmasker/hard
cp $wkdir/$GENOME $wkdir/${GENOME%%.*}/repeatmasker/hard/
RepeatMasker -s -lib $wkdir/${GENOME%%.*}/blastx/${GENOME%%.*}.clust.class.noprot.fa -gff -pa 16 $wkdir/${GENOME%%.*}/repeatmasker/hard/$GENOME
