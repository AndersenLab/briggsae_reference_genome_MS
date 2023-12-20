#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="tpsi"

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"
tpdir="/projects/b1059/projects/Nicolas/software/tPSI_deploy"

source activate tPSI

cd $wkdir/${file%%.*}/transposonPSI
perl $tpdir/transposonPSI.pl $wkdir/$file nuc
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,$9}' ${file%.*}.fa.TPSI.allHits.chains.bestPerLocus.gff3 > ${file%.*}.TPSI.allHits.chains.bestPerLocus.bed
bedtools getfasta -name -fi $wkdir/$file -bed ${file%.*}.TPSI.allHits.chains.bestPerLocus.bed -fo ${file%.*}.TPSI.allHits.chains.bestPerLocus.fa 
awk '!/^>/ { next } { getline seq } length(seq) >= 50 { print $0 "\n" seq }' ${file%.*}.TPSI.allHits.chains.bestPerLocus.fa  > ${file%.*}.TPSI.allHits.chains.bestPerLocus_50bp.fa 

