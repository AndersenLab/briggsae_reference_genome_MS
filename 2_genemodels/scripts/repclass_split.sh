#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="repclass_s"

source activate rep_mod

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"


cp $wkdir/${file%%.*}/vsearch/${file%%.*}.repeats.unclass.clust.fa $wkdir/${file%%.*}/classify
cd $wkdir/${file%%.*}/classify
mv ${file%%.*}.repeats.unclass.clust.fa partitioned.repeats.unclass.clust.fa
perl $wkdir/scripts/fasta-splitter.pl --n-parts 250 --line-length 0 --measure count $wkdir/${file%%.*}/classify/partitioned.repeats.unclass.clust.fa 
rm partitioned.repeats.unclass.clust.fa
for i in *.fa; do RepeatClassifier -consensi $i -engine ncbi; done
cat *.fa.classified > ${file%%.*}.clust.class.fa


sbatch --export=CLASSREP=${file%%.*}.clust.class.fa --output=$wkdir/scripts/logs/${file%%.*}.blastx.oe $wkdir/scripts/blastx.sh
