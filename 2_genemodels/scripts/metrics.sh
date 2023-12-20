#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --output="metrics.oe"
#SBATCH --job-name="metrics"


source activate customR

wkdir="/projects/b1042/AndersenLab/work/nic/N2"

#Rscript aminoacids elegansprot blastout

cd $wkdir/gff2prot/braker/
Rscript $wkdir/scripts/metrics.R ${GENOME%%.*}.braker.aa.fa $wkdir/elegansprot/elegans.protein.clean.fa ${GENOME%%.*}.braker.pb.out

cd $wkdir/gff2prot/genemark/
Rscript $wkdir/scripts/metrics.R ${GENOME%%.*}.genemark.aa.fa $wkdir/elegansprot/elegans.protein.clean.fa ${GENOME%%.*}.genemark.pb.out

cd $wkdir/gff2prot/maker/
Rscript $wkdir/scripts/metrics.R ${GENOME%%.*}.maker.aa.fa $wkdir/elegansprot/elegans.protein.clean.fa ${GENOME%%.*}.maker.pb.out

cd $wkdir/gff2prot/augustus/nohints/
Rscript $wkdir/scripts/metrics.R ${GENOME%%.*}.augustus.nhints.aa.fa $wkdir/elegansprot/elegans.protein.clean.fa ${GENOME%%.*}.augustus.nhints.pb.out

cd $wkdir/gff2prot/augustus/default/
Rscript $wkdir/scripts/metrics.R ${GENOME%%.*}.augustus.whints.aa.fa $wkdir/elegansprot/elegans.protein.clean.fa ${GENOME%%.*}.augustus.whints.pb.out

