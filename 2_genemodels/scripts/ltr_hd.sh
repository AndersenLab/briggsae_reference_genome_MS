#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="ltrharvest"


wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"
source activate ltr_hd

cd $wkdir/${file%%.*}/ltr_hd
gt suffixerator -db $wkdir/$file -indexname ${file%%.*} -suf -lcp -des -ssp -sds -dna -lossless -v
gt ltrharvest -index ${file%%.*} -seqids yes -md5 yes -tabout no | gt gff3 -sort -retainids -tidy > ${file%%.*}.ltrharvest.gff3
gt ltrharvest -index ${file%%.*} -out ${file%%.*}.ltr.fa
gt ltrdigest -hmms $wkdir/scripts/all_hmms/*.hmm -encseq ${file%%.*} -v < ${file%%.*}.ltrharvest.gff3 > ${file%%.*}.ltrdigest.gff3
gt select -rule_files $wkdir/scripts/filter_protein_match.lua -- < ${file%%.*}.ltrdigest.gff3 > filtered.${file%%.*}.ltrdigest.gff3
gt extractfeat -type LTR_retrotransposon -encseq ${file%%.*} filtered.${file%%.*}.ltrdigest.gff3 > filtered.${file%%.*}.ltrdigest.fa

