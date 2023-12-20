#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="repmodeler"

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"

source activate rep_mod

echo $GENOME
cd $wkdir/${GENOME%%.*}/repeatmodeler
BuildDatabase -name ${GENOME%.*} $wkdir/$GENOME 
RepeatModeler -engine ncbi -pa 24 -database ${GENOME%.*} > ${GENOME%.*}.rm.out
