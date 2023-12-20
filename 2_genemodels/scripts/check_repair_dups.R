args <- commandArgs(trailingOnly = TRUE)

library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(zoo)

#read gff
file <- args[1]
gff <- read.gff(file)

#extract all L1s
L1_features <- gff %>% dplyr::filter(type=="gene") 

#extract L2s and separate features
L2_features <- gff %>% dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,c("Feature_ID","Parent_ID","Other"),sep=";",extra="merge") 
L2_features$Parent_ID <- gsub(".*=","ID=",L2_features$Parent_ID)

#L2_check <- L2_features %>% dplyr::filter(Parent_ID %in% L1_features_pc$Feature_ID) #<- this was used to check that mRNA parents are protein-coding genes (duh!)

#extract L3 features
L3_features <- gff %>% dplyr::filter(grepl("exon|intron|CDS|start_codon|stop_codon|five_prime_UTR|three_prime_UTR",type)) %>%
  tidyr::separate(attributes,c("Feature_ID","Parent_ID","Other"),sep=";",extra="merge")
L3_features$Parent_ID <- gsub(".*=","ID=",L3_features$Parent_ID)

#exctract and summarise length of CDS features per transcript
L2_CDS_lengths <- L3_features %>% dplyr::filter(type=="CDS") %>%
  dplyr::mutate(cds_length=end-start) %>%
  dplyr::group_by(Parent_ID) %>%
  # dplyr::mutate(ORF_length=sum(cds_length)) %>%
  dplyr::summarise(CDSlengths=toString(cds_length)) %>%
  dplyr::ungroup() 

#extract and summarise start positions of every CDS in each transcript
L2_CDS_starts <- L3_features %>% dplyr::filter(type=="CDS") %>%
  dplyr::group_by(Parent_ID) %>%
  dplyr::summarise(CDSstarts=toString(start)) %>%
  dplyr::ungroup()

#extract and summarise length of intron features per transcript
L2_introns_lengths <- L3_features %>% dplyr::filter(type=="exon") %>%
  dplyr::arrange(Parent_ID,seqid,start)%>%
  dplyr::group_by(Parent_ID) %>%
  dplyr::mutate(intron=lead(start)-end) %>%
  dplyr::select(-Other,-score,-phase) %>%
  dplyr::ungroup()%>%
  replace(is.na(.), -1) %>%
  dplyr::filter(intron!=-1) %>%
  dplyr::group_by(Parent_ID) %>%
  dplyr::summarise(intron_lengths=toString(intron))

#extact and summarise start positions of every intron in each transcript
L2_intron_starts <- L3_features %>% dplyr::filter(type=="exon") %>%
  dplyr::arrange(Parent_ID,seqid,start)%>%
  dplyr::group_by(Parent_ID) %>%
  dplyr::mutate(intron_start=lag(end)) %>%
  dplyr::select(-Other,-score,-phase) %>%
  dplyr::ungroup()%>%
  replace(is.na(.), -1) %>%
  dplyr::filter(intron_start!=-1) %>%
  dplyr::group_by(Parent_ID) %>%
  dplyr::summarise(intron_starts=toString(intron_start))
 
#append L3 feature (intron, CDS) stats extracted above to L2 features (mRNA)
L2_L3stats <- dplyr::left_join(L2_CDS_lengths,L2_CDS_starts,by="Parent_ID")  %>% 
  dplyr::left_join(., L2_introns_lengths, by='Parent_ID') %>%
  dplyr::left_join(., L2_intron_starts, by='Parent_ID') 

colnames(L2_L3stats) <- c("Feature_ID","CDS_lengths","CDS_starts","intron_lengths","intron_starts")

#identify BRAKER trascripts with redundant CDS when compared with StringTie/Transdecoder transcripts in the same locus
L2_features_wORFlen <- dplyr::left_join(L2_features,L2_L3stats) %>%
  dplyr::mutate(source2=ifelse(source=="transdecoder","TRANSDECODER","BRAKER")) %>%
  dplyr::group_by(Parent_ID,CDS_lengths,CDS_starts) %>%
  dplyr::mutate(group_size=n()) %>%
  dplyr::mutate(element_ID=1:n()) %>%
  dplyr::mutate(class=ifelse(all(source2=="BRAKER"),"ALL_BRAKER",ifelse(all(source2=="TRANSDECODER"),"ALL_TD","MIX"))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(group_size>1) %>%
  dplyr::mutate(tran_len=end-start) %>%
  dplyr::arrange(Parent_ID,desc(source),desc(tran_len)) 

BRAKER_toRM <- L2_features_wORFlen %>% dplyr::filter(class=="MIX") %>% dplyr::filter(source2=="BRAKER")

#check if any StringTie/Transdecoder transcripts with CDS redundancies have intron chain redundancies 
remain_TD <- L2_features_wORFlen %>% dplyr::filter(class=="MIX") %>% 
  dplyr::filter(source2=="TRANSDECODER") %>%
  dplyr::group_by(Parent_ID,CDS_lengths,CDS_starts) %>%
  dplyr::mutate(new_group_size=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(new_group_size>1) %>%
  dplyr::select(-new_group_size,-class)

#It appears that there are no redundant StringTie transcripts - this dataframe is not used
all_TD <- dplyr::bind_rows(L2_features_wORFlen %>% dplyr::filter(class=="ALL_TD") %>% dplyr::select(-class),remain_TD) %>%
  dplyr::group_by(Parent_ID,CDS_lengths,CDS_starts,intron_lengths,intron_starts) %>%
  dplyr::mutate(new_group_size=n()) %>%
  dplyr::ungroup() 

#get list of of genes affected by redundancies
affectedL1s <- BRAKER_toRM %>% dplyr::distinct(Parent_ID)

#filter out redundant L2 features
remaining_L2_features <- L2_features %>% dplyr::filter(!Feature_ID %in% BRAKER_toRM$Feature_ID)
#filter out redundant L3 features
reamining_L3_features <- L3_features %>% dplyr::filter(!Parent_ID %in% BRAKER_toRM$Feature_ID)

#reset GFF format to standard
reamining_L3_features$Parent_ID <- gsub(".*=","Parent=",reamining_L3_features$Parent_ID)
reamining_L3_features_cleanGFF <- reamining_L3_features %>% tidyr::unite("attributes",Feature_ID,Parent_ID,Other,sep=";",na.rm = T)
remaining_L2_features$Parent_ID <- gsub(".*=","Parent=",remaining_L2_features$Parent_ID)
remaining_L2_features_cleanGFF <- remaining_L2_features %>% tidyr::unite("attributes",Feature_ID,Parent_ID,Other,sep=";",na.rm = T)

#combine all filtered features into final GFF
clean_GFF <- bind_rows(L1_features,remaining_L2_features_cleanGFF,reamining_L3_features_cleanGFF) %>% dplyr::arrange(seqid,start)

#get dir and filename
dir <- dirname(file)
fileName <- basename(file)

#set output filename
out1 <- gsub(".gff",".dedup.gff",fileName)
outName1 <- paste0(dir,"/",out1)

#write GFF
write.table(clean_GFF,file=outName1,quote = F,sep = "\t",row.names = F,col.names = F,na = ".")

