#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

#SBATCH --job-name="blastx"

source activate blastx

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"

cd $wkdir/${CLASSREP%%.*}/blastx
blastx -db $wkdir/wormbase_dbs/elegansprot/elegans.protein -num_threads 24 -outfmt 6 -query $wkdir/${CLASSREP%%.*}/classify/$CLASSREP -out BLAST_${CLASSREP%%.*}.res.txt -evalue 0.001

sbatch --export=file=${CLASSREP%%.*}.genome.fa --output=$wkdir/scripts/logs/${CLASSREP%%.*}.curate.oe $wkdir/scripts/curate.sh
