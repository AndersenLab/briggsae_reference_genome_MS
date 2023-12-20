### Author: Nicolas Moya - Andersen Lab
### The aim of this script is to merge BRAKER and StringTie/Transdecoder gene predictions into a single, non-overlaping, non-redundant gene set

#working directory
$wkdir=/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions

cd $wkdir/QX1410/predictions/merger/masked_merger/final_merger


# merge BRAKER and StringTie models
agat_sp_merge_annotations.pl -f $wkdir/QX1410/predictions/braker/braker2.1.6_masked/braker.gff3 -f $wkdir/QX1410/predictions/merger/masked_merger/stringtie_filtered/QX1410.stringtie.clean_merger.filtered.final_kept.gff -o QX1410.SB.merger1.gff

#fix overlappping genes with shared CDS
agat_sp_fix_overlaping_genes.pl -f QX1410.SB.merger1.gff -o QX1410.SB.merger1.fixed_CDSoverlaps.gff

#identify non-CDS overlaps - primary output is the coordinates of overlaps (parents.coordinates.txt)
Rscript find_overlap_coordinates.R QX1410.SB.merger1.fixed_CDSoverlaps.gff

#filter out genes that fall within overlap coordinates in merger GFF
agat_sp_filter_record_by_coordinates.pl -gff QX1410.SB.merger1.fixed_CDSoverlap.gff -c QX1410.SB.merger1.fixed_CDSoverlap.parents.coordinates.txt -o QX1410.SB.merger.overlaps

#apply same filter to BRAKER GFF - we will replace the overlaps in merger GFF with original BRAKER models
agat_sp_filter_record_by_coordinates.pl -gff $wkdir/QX1410/predictions/braker/braker2.1.6_masked/braker.gff -c QX1410.merger.fixedCDS.parents.coordinates.txt -o QX1410.braker.overlaps

#concatenate all files that contain features within coordinates in BRAKER GFF
cd QX1410.braker.overlaps
mkdir within_coords
mv *.gff3 within_coords/
mv within_coords/remaining.gff3 $PWD
cd within_coords
#concatenate files
cat *.gff3 | grep -v "^#" > QX1410.braker.within_coordinates.gff
#fix pre-existing CDS overlaps in BRAKER GFF 
agat_sp_fix_overlaping_genes.pl -f QX1410.braker.within_coordinates.gff -o QX1410.braker.within_coordinates.fixedCDS_overlap.gff 
cd ../../

#organize files within merger GFF filter ouput for consistency
cd QX1410.SB.merger.overlaps
mkdir within_coords
mv *.gff3 within_coords/
mv within_coords/remaining.gff3 $PWD
cd ../


#merge overlap_filtered
agat_sp_merge_annotations.pl -f QX1410.SB.merger.overlaps/remaining.gff3 -f QX1410.braker.overlaps/within_coords/QX1410.braker.within_coordinates.fixedCDS_overlap.gff -o QX1410.SB.replaced_coords.gff

#check for any remaining non-CDS overlaps (carried over from pre-existing non-CDS overlaps in BRAKER GFF)
Rscript find_overlap_coordinates.R QX1410.SB.replaced_coords.gff

#get overlapped transcript list to remove
grep -v "seqid" QX1410.SB.replaced_coords.donors.txt | awk '{print $10}' | sed 's/ID=//' > QX1410.genes_toRM.txt

#'nbis' named genes have different Parent/Child structure, so we have to treat them differently in order to remove them properly from the gff
grep "nbis" QX1410.genes_toRM.txt > QX1410.nbis_genes_toRM.txt
grep -wFf QX1410.nbis_genes_toRM.txt QX1410.SB.replaced_coords.gff | grep "mRNA" | awk '{print $9}' | sed 's/;.*//' | sed 's/ID=//' > QX1410.nbis_tran_toRM.txt
grep -v "nbis" QX1410.genes_toRM.txt > QX1410.normal_genes_toRM.txt
cat QX1410.nbis_genes_toRM.txt QX1410.nbis_tran_toRM.txt QX1410.normal_genes_toRM.txt > QX1410.genes_toRM.clean.txt

#remove overlaps 
grep -vwFf QX1410.genes_toRM.clean.txt QX1410.SB.replaced_coords.gff > QX1410.SB.final_merger.noOverlaps.clean.gff

#rename features IDs with cohesive prefixes

agat_sp_manage_IDs.pl -gff QX1410.SB.final_merger.noOverlaps.clean.gff --prefix QX1410. --type_dependent --tair -o QX1410.SB.final_merger.renamed.gff

#remove transcripts with redundant CDS and intron chain (output will be .dedup.gff)
Rscript check_repair_dups.R QX1410.SB.final_merger.renamed.gff

#since we removed transcripts, we'll rename once again. A slightly different prefix will be used due to a limitation in AGAT script.

agat_sp_manage_IDs.pl -gff QX1410.SB.final_merger.renamed.dedup.gff --prefix QX1410.g --type_dependent --tair -o QX1410.SB.final_anno.dedup.renamed.gff

#final annotation was copied to '~/final_annotations/' directory in the GitHub repository.

