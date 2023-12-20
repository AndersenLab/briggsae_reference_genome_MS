### Author: Nicolas Moya - Andersen Lab
### The aim of this script is to do reciprocal BLAST between C. elegans N2 and C. briggsae (QX1410, VX34, AF16) protein FASTAs.
### libraries were generated using 'makeblastdb' from the BLAST suite, using BLAST version 5 and dbtype 'prot'

wkdir="/projects/b1059/projects/Nicolas/c.briggsae/blast_analysis"

query=$1
library=$2
libname="$(basename $library)"

#BLAST C.b. protein FASTA against C.e. N2 library (forward)
cd $wkdir/blast_out/forward
blastp -query $query -db $wkdir/N2_library/c_elegans.PRJNA13758.WS279.protein_coding.prot.fa -out ${query%.*}.pb.out -outfmt 6 -num_threads 24

cd $wkdir/blast_out/reciprocal
#BLAST C.e. N2 protein FASTA against C.b. library (reciprocal)
blastp -query $wkdir/N2_library/c_elegans.PRJNA13758.WS279.protein_coding.prot.fa -db $library -out ${libname%.*}.recipro.pb.out -outfmt 6 -num_threads 24