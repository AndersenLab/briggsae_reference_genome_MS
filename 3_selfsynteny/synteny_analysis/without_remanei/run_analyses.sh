# fix GFFs for CINOP, CNIGO, CREMA
echo "Fixing GFFs..."
bash ../scripts/fix_CINOP.sh gffs/c_inopinata.PRJDB5687.WS279.annotations.gff3_longest_isoforms
bash ../scripts/fix_CNIGO.sh gffs/c_nigoni.PRJNA384657.WS279.annotations.gff3_longest_isoforms

# print TSV
echo "Printing TSVs..."
python3 ../scripts/print_tsv.py -o Orthogroups.single_copy.txt -1 CELEG -2 CBRIG -y gffs/c_elegans.PRJNA13758.WS279.annotations.gff3_longest_isoforms -z gffs/caenorhabditis_briggsae_QX1410_v2.1_annotations.nuclear.gff3 | sort -g -k 6 >CELEG_CBRIG.txt
python3 ../scripts/print_tsv.py -o Orthogroups.single_copy.txt -1 CINOP -2 CNIGO -y gffs/c_inopinata.PRJDB5687.WS279.annotations.gff3_longest_isoforms -z gffs/c_nigoni.PRJNA384657.WS279.annotations.gff3_longest_isoforms | sort -g -k 6 >CINOP_CNIGO.txt

# calculate synteny
echo "Calculating synteny..."
python3 ../scripts/calculate_synteny.py -t CELEG_CBRIG.txt | grep -v Total >CELEG_CBRIG.syn
python3 ../scripts/calculate_synteny.py -t CINOP_CNIGO.txt | grep -v Total >CINOP_CNIGO.syn

# make syteny by chromosome plot 
Rscript ../scripts/plot_synteny_by_chr.R CELEG_CBRIG.syn CINOP_CNIGO.syn

# calculate synteny windows 
python3 ../scripts/synteny_windows.py -t CELEG_CBRIG.txt -w 500000 >CELEG_CBRIG.syn_windows
python3 ../scripts/synteny_windows.py -t CINOP_CNIGO.txt -w 500000 >CINOP_CNIGO.syn_windows

# identify synteny blocks 
python3 ../scripts/synteny_blocks.py CELEG_CBRIG.txt >CELEG_CBRIG.syn_blocks
python3 ../scripts/synteny_blocks.py CINOP_CNIGO.txt >CINOP_CNIGO.syn_blocks

# make windows + rearrangement plots 
Rscript ../scripts/plot_syteny_windows_rearrangements.R
