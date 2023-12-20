### Author: Nicolas Moya - Andersen Lab
### The aim of this script is to filter out incomplete and abnormal transcripts from StringTie/TransDecoder gene predictions.

#working directory
$wkdir=/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions


#reformat FASTA to be bioperl compliant for use with AGAT (60 characters per line)
#fasta formatter tool from FASTX-Toolkit (conda install -c bioconda fastx_toolkit)
fasta_formatter -i $wkdir/genomes/QX1410/decompressed/QX1410.genome.v2.1.fa -w 60 > $wkdir/genomes/QX1410/decompressed/QX1410.genome.v2.1.rf.fa


cd $wkdir/QX1410/predictions/merger/masked_merger/stringtie_filtered

#agat_sp_filter_incomplete_gene_coding_models.pl -gff $wkdir/$strain/predictions/stringtie/$strain.stringtie.TDm70.gff3 --fasta $wkdir/genomes/$strain/decompressed/$strain.genome.v2.1.rf.fa -o $strain.stringtie.CDSfiltered.gff3
agat_sp_filter_incomplete_gene_coding_models.pl -gff cd $wkdir/QX1410/predictions/stringtie/QX1410.stringtie.TDm70.gff3 --fasta cd $wkdir/genomes/QX1410/decompressed/QX1410.genome.v2.1.rf.fa --ct 0 -o QX1410.stringtie.filtered.gff3


#fix CDS of incomplete models (codon table 0 to force only AUG start codons)
agat_sp_fix_longestORF.pl -gff QX1410.stringtie.filtered_incomplete.gff3  -f $wkdir/genomes/QX1410/decompressed/QX1410.genome.v2.1.rf.fa --ct 0 -m 1,3,4,5 -o QX1410.stringtie.filtered_incomplete.fixedCDS.gff

#manually remove problematic transcript(s) with no CDS - this issue was spotted downstream when filtering incomplete models
grep -vwFf manual_removals.txt QX1410.stringtie.filtered_incomplete.fixedCDS-only_modified.gff  > QX1410.stringtie.filtered_incomplete.fixedCDS-only_modified2.gff 

#merge with complete stringtie models
agat_sp_merge_annotations.pl -f QX1410.stringtie.filtered.gff -f QX1410.stringtie.filtered_incomplete.fixedCDS-only_modified2.gff -o QX1410.stringtie.clean_merger.gff 

#remove transcripts that are incomplete after CDS fix
agat_sp_filter_incomplete_gene_coding_models.pl -gff QX1410.stringtie.clean_merger.gff -fa $wkdir/genomes/QX1410/decompressed/QX1410.genome.v2.1.rf.fa --ct 0 -o QX1410.stringtie.clean_merger.filtered.gff

#remove abnormal transcripts (>1 non-coding exon)
Rscript filter_abnormal_transcripts.R QX1410.stringtie.clean_merger.filtered.gff 0

#final StringTie/TransDecoder models will be in 'QX1410.stringtie.clean_merger.filtered.final_kept.gff'
