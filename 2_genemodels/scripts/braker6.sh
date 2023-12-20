#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="braker6"
#SBATCH --output="braker6.oe"

source activate braker6_env
module load samtools

/projects/b1059/projects/Nicolas/software/braker-2.1.6/BRAKER/scripts/braker.pl --gff3 --useexisting --species=QX1410.b6.v2.1.masked --genome=/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/QX1410/repeatmasker/soft/QX1410.genome.v2.1.fa.masked --bam=/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/QX1410/alignments/Aligned.sortedByCoord.out.bam --AUGUSTUS_SCRIPTS_PATH=/home/ndm7437/.conda/envs/braker6_env/bin/ --cores=24

