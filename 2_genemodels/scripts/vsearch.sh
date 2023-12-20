#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="vsearch"

source activate ltr_hd

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"


cp $wkdir/${file%%.*}/ltr_hd/filtered.${file%%.*}.ltrdigest.fa $wkdir/${file%%.*}/vsearch
cp $wkdir/${file%%.*}/transposonPSI/${file%.*}.TPSI.allHits.chains.bestPerLocus_50bp.fa $wkdir/${file%%.*}/vsearch
cp $wkdir/${file%%.*}/repeatmodeler/RM_*/consensi.fa $wkdir/${file%%.*}/vsearch
cp $wkdir/scripts/rhabrep/rhabditida_repeats.fa $wkdir/${file%%.*}/vsearch
cp $wkdir/scripts/rhabrep/dfam.repeats.fa $wkdir/${file%%.*}/vsearch
cd $wkdir/${file%%.*}/vsearch
cat * > ${file%%.*}.repeats.unclass.fa
vsearch --cluster_fast ${file%%.*}.repeats.unclass.fa -consout ${file%%.*}.repeats.unclass.clust.fa -msaout ${file%%.*}.repeats.unclass.aligned.fa --id 0.8

sbatch --export=file=$file --output=$wkdir/scripts/logs/${file%%.*}.repclass.oe $wkdir/scripts/repclass_split.sh
