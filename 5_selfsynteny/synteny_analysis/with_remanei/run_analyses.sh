# fix GFFs for CINOP, CNIGO, CREMA
echo "Fixing GFFs..."
#bash ../scripts/fix_CINOP.sh gffs/c_inopinata.PRJDB5687.WS279.annotations.gff3_longest_isoforms
#bash ../scripts/fix_CNIGO.sh gffs/c_nigoni.PRJNA384657.WS279.annotations.gff3_longest_isoforms
#bash ../scripts/fix_CREMA.sh gffs/GCA_010183535.1_CRPX506_genomic.gff3_longest_isoforms

# print TSV
echo "Printing TSVs..."
#python3 ../scripts/print_tsv.py -o Orthogroups.single_copy.txt -1 CBRIG -2 CREMA -y gffs/caenorhabditis_briggsae_QX1410_v2.1_annotations.nuclear.gff3 -z gffs/GCA_010183535.1_CRPX506_genomic.gff3_longest_isoforms | sort -g -k 6 >CBRIG_CREMA.txt
#python3 ../scripts/print_tsv.py -o Orthogroups.single_copy.txt -1 CNIGO -2 CREMA -y gffs/c_nigoni.PRJNA384657.WS279.annotations.gff3_longest_isoforms -z gffs/GCA_010183535.1_CRPX506_genomic.gff3_longest_isoforms | sort -g -k 6 >CNIGO_CREMA.txt
#python3 ../scripts/print_tsv.py -o Orthogroups.single_copy.txt -1 CELEG -2 CREMA -y gffs/c_elegans.PRJNA13758.WS279.annotations.gff3_longest_isoforms -z gffs/GCA_010183535.1_CRPX506_genomic.gff3_longest_isoforms | sort -g -k 6 >CELEG_CREMA.txt
#python3 ../scripts/print_tsv.py -o Orthogroups.single_copy.txt -1 CINOP -2 CREMA -y gffs/c_inopinata.PRJDB5687.WS279.annotations.gff3_longest_isoforms -z gffs/GCA_010183535.1_CRPX506_genomic.gff3_longest_isoforms | sort -g -k 6 >CINOP_CREMA.txt

# calculate synteny
echo "Calculating synteny..."
python3 ../scripts/calculate_synteny.py -t CBRIG_CREMA.txt | grep -v Total >CBRIG_CREMA.syn
python3 ../scripts/calculate_synteny.py -t CNIGO_CREMA.txt | grep -v Total | grep -v PDUG >CNIGO_CREMA.syn
python3 ../scripts/calculate_synteny.py -t CELEG_CREMA.txt | grep -v Total >CELEG_CREMA.syn
python3 ../scripts/calculate_synteny.py -t CINOP_CREMA.txt | grep -v Total >CINOP_CREMA.syn

# plot
echo "Making plots..."
Rscript ../scripts/plot_synteny_by_chr.R CBRIG_CREMA.syn CNIGO_CREMA.syn CBRIG_CNIGO_synteny_by_chr.pdf
Rscript ../scripts/plot_synteny_by_chr.R CELEG_CREMA.syn CINOP_CREMA.syn CELEG_CINOP_synteny_by_chr.pdf
Rscript ../scripts/comparison_with_remanei.R
