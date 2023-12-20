#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --output=busco.oe
#SBATCH --job-name="busco"

source activate busco
wkdir="/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions"

busco -i $file -c 12 -m prot -l /projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/scripts/busco_downloads/lineages/nematoda_odb10 -o ${file%.*}.busco
