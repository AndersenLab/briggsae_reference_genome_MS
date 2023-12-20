#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="transdec"

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"

source activate tama

gtf_genome_to_cdna_fasta.pl $file $wkdir/${file%%.*}.genome.v2.1.fa > ${file%%.*}.transcripts.fa
gtf_to_alignment_gff3.pl $file > ${file%%.*}.transcripts.gff3
TransDecoder.LongOrfs -t ${file%%.*}.transcripts.fa -m 70
TransDecoder.Predict -t ${file%%.*}.transcripts.fa

cdna_alignment_orf_to_genome_orf.pl \
     ${file%%.*}.transcripts.fa.transdecoder.gff3 \
     ${file%%.*}.transcripts.gff3 \
     ${file%%.*}.transcripts.fa > ${file%%.*}.stringtie.TDm70.gff3
