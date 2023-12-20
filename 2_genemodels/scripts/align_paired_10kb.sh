#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="STARalign"

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"
source activate star

cd $wkdir/${GENOME%%.*}/alignments/
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 600000000000 \
--genomeDir . \
--genomeFastaFiles $wkdir/$GENOME \
--genomeSAindexNbases 12 \
--alignIntronMax 10000
STAR \
--runThreadN 24 \
--genomeDir . \
--outSAMtype BAM Unsorted SortedByCoordinate \
--twopassMode Basic \
--readFilesCommand zcat \
--alignIntronMax 10000 \
--readFilesIn $wkdir/shortreads/${GENOME%%.*}/${GENOME%%.*}.reads.f.fq.gz $wkdir/shortreads/${GENOME%%.*}/${GENOME%%.*}.reads.r.fq.gz

