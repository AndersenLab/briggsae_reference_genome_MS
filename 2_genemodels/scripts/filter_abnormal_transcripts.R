args <- commandArgs(trailingOnly = TRUE)

library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(zoo)

### pull CDS2exon ratios for a given gff
pull_ratios <-function(gff,isWB) {
  
  if(isWB==0) {
  stringtie_L1_features <- gff %>% dplyr::filter(type=="gene") %>% 
    tidyr::separate(attributes,into=c("Feature_ID","Name"),sep=";")
  
  stringtie_L2_features <- gff %>% dplyr::filter(type=="mRNA")%>% 
    tidyr::separate(attributes,into=c("Feature_ID","Parent_ID","Name"),sep=";")
  stringtie_L2_features$Parent_ID <- gsub(".*=","ID=",stringtie_L2_features$Parent_ID)
  
  stringtie_L3_features <- gff %>% dplyr::filter(grepl("exon|intron|CDS|start_codon|stop_codon|five_prime_UTR|three_prime_UTR",type))%>% 
    tidyr::separate(attributes,into=c("Feature_ID","Parent_ID"),sep=";")
  stringtie_L3_features$Parent_ID <- gsub(".*=","ID=",stringtie_L3_features$Parent_ID)
  
  L1L2_relationships <- stringtie_L2_features %>% dplyr::select(Feature_ID,Parent_ID)
  colnames(L1L2_relationships) <- c("mRNA_ID","Gene_ID")
  
  stringtie_L3_features_wL1 <- dplyr::left_join(stringtie_L3_features,L1L2_relationships,by=c("Parent_ID"="mRNA_ID"))
  
  } else {
    stringtie_L1_features <- gff %>% dplyr::filter(type=="gene") %>% 
      tidyr::separate(attributes,into=c("Feature_ID","Name"),sep=";")
    
    stringtie_L2_features <- gff %>% dplyr::filter(type=="mRNA")%>% 
      tidyr::separate(attributes,into=c("Feature_ID","Parent_ID"),sep=";")
    stringtie_L2_features$Parent_ID <- gsub(".*=","ID=",stringtie_L2_features$Parent_ID)
    
    stringtie_L3_features_nonCDS <- gff %>% dplyr::filter(grepl("exon|intron|start_codon|stop_codon|five_prime_UTR|three_prime_UTR",type)) %>%
      tidyr::separate(attributes,c("Parent_ID","Other"),sep=";",extra="merge") %>% 
      dplyr::select(-Other)
    stringtie_L3_features_nonCDS$Parent_ID <- gsub(".*=","ID=",stringtie_L3_features_nonCDS$Parent_ID) 
      
    
    stringtie_CDS <- gff %>% dplyr::filter(type=="CDS") %>% 
      tidyr::separate(attributes,c("Feature_ID","Parent_ID","Other"),sep=";",extra="merge") %>% 
      dplyr::select(-Other,-Feature_ID)
    
    multiTranscript_CDS <- stringtie_CDS %>% dplyr::filter(grepl(",",Parent_ID)) %>%
      tidyr::separate_rows(Parent_ID,sep = ",")
    
    missingParent <- multiTranscript_CDS %>% dplyr::filter(!grepl("Parent",Parent_ID))
    missingParent$Parent_ID <- gsub("Transcript","Parent=Transcript",missingParent$Parent_ID) 
    
    fixed_CDS <- bind_rows(missingParent,multiTranscript_CDS %>% dplyr::filter(grepl("Parent",Parent_ID)),stringtie_CDS %>% dplyr::filter(!grepl(",",Parent_ID)))
    fixed_CDS$Parent_ID <- gsub(".*=","ID=",fixed_CDS$Parent_ID)
    fixed_L3s <- bind_rows(stringtie_L3_features_nonCDS,fixed_CDS) %>% dplyr::arrange(seqid,start)
    
    #L3_check <- dplyr::filter(fixed_L3s,Parent_ID %in% stringtie_L2_features$Feature_ID)
    
    L1L2_relationships <- stringtie_L2_features %>% dplyr::select(Feature_ID,Parent_ID)
    colnames(L1L2_relationships) <- c("mRNA_ID","Gene_ID")
    
    stringtie_L3_features_wL1 <- dplyr::left_join(fixed_L3s,L1L2_relationships,by=c("Parent_ID"="mRNA_ID"))
  }
  
  stringtie_L2_exon_lengths <- stringtie_L3_features_wL1 %>% dplyr::filter(type=="exon") %>%
    dplyr::mutate(exon_length=end-start) %>%
    dplyr::group_by(Parent_ID) %>%
    dplyr::mutate(exon_count=n()) %>%
    dplyr::mutate(exon_sum=sum(exon_length)) %>%
    dplyr::distinct(Parent_ID,.keep_all = T) %>% 
    dplyr::ungroup() %>%
    dplyr::select(Parent_ID,exon_sum,exon_count)
  
  stringtie_L2_UTR_lengths <- stringtie_L3_features_wL1 %>% dplyr::filter(grepl("five_prime_UTR|three_prime_UTR",type)) %>%
    dplyr::mutate(UTR_length=end-start) %>%
    dplyr::group_by(Parent_ID) %>%
    dplyr::mutate(three_prime_count=sum(type=="three_prime_UTR")) %>%
    dplyr::mutate(five_prime_count=sum(type=="five_prime_UTR")) %>%
    dplyr::mutate(UTR_sum=sum(UTR_length)) %>%
    dplyr::distinct(Parent_ID,.keep_all = T) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(normal_UTR=ifelse(five_prime_count+three_prime_count < 3 ,"Y","N")) %>%
    dplyr::mutate(UTR_counts=five_prime_count+three_prime_count) %>%
    dplyr::select(Parent_ID,normal_UTR,UTR_counts,UTR_sum)
  
  stringtie_L2_CDS_lengths <- stringtie_L3_features_wL1 %>% dplyr::filter(type=="CDS") %>%
    dplyr::mutate(CDS_length=end-start) %>%
    dplyr::group_by(Parent_ID) %>%
    dplyr::mutate(CDS_sum=sum(CDS_length)) %>%
    dplyr::mutate(CDS_count=n()) %>%
    dplyr::distinct(Parent_ID,.keep_all = T) %>% 
    dplyr::ungroup() %>%
    dplyr::select(Parent_ID,CDS_sum,CDS_count)
  
  #check <- stringtie_L3_features_wL1 %>% group_by(Parent_ID) %>% dplyr::mutate(count=sum(type=="three_prime_UTR" | type=="five_prime_UTR"))
  
  stringtie_all_ratio <- dplyr::left_join(stringtie_L2_exon_lengths,stringtie_L2_CDS_lengths,by="Parent_ID") %>% dplyr::left_join(.,stringtie_L2_UTR_lengths,by="Parent_ID") %>% dplyr::mutate(ratio=CDS_sum/exon_sum) 
  #colnames(stringtie_all_ratio) <- c("Gene_ID","Transcript_ID","ratio")
  
  stringtie_UTR2exon_normalUTR_ratio <-    stringtie_all_ratio %>% dplyr::filter(normal_UTR=="Y")
  stringtie_UTR2exon_abnormalUTR_3c_ratio <-   stringtie_all_ratio %>% dplyr::filter(normal_UTR=="N" & UTR_counts==3)
  stringtie_UTR2exon_abnormalUTR_4c_ratio <-   stringtie_all_ratio %>% dplyr::filter(normal_UTR=="N" & UTR_counts==4)
  stringtie_UTR2exon_abnormalUTR_5c_ratio <-   stringtie_all_ratio %>% dplyr::filter(normal_UTR=="N" & UTR_counts==5)
  stringtie_UTR2exon_abnormalUTR_6p_ratio <-   stringtie_all_ratio %>% dplyr::filter(normal_UTR=="N" & UTR_counts>5)
  
  all_ratios <- bind_rows(stringtie_all_ratio %>% dplyr::select(ratio),
                          stringtie_UTR2exon_normalUTR_ratio%>%dplyr::select(ratio),
                          stringtie_UTR2exon_abnormalUTR_3c_ratio%>%dplyr::select(ratio),
                          stringtie_UTR2exon_abnormalUTR_4c_ratio%>%dplyr::select(ratio),
                          stringtie_UTR2exon_abnormalUTR_5c_ratio%>%dplyr::select(ratio),
                          stringtie_UTR2exon_abnormalUTR_6p_ratio%>%dplyr::select(ratio),
                          .id="source")
  
  all_ratios <- dplyr::mutate(all_ratios,source2=case_when(source==1 ~ "All Transcripts",
                                                           source==2 ~ "Normal Transcripts (1-2 UTRs)",
                                                           source==3 ~ "Abnormal Transcripts (3 UTRs)",
                                                           source==4 ~ "Abnormal Transcripts (4 UTRs)",
                                                           source==5 ~ "Abnormal Transcripts (5 UTRs)",
                                                           source==6 ~ "Abnormal Transcripts (6+ UTRs)"))
  all_ratios$source2 <- factor(all_ratios$source2, levels = c("All Transcripts",
                                                            "Normal Transcripts (1-2 UTRs)",
                                                            "Abnormal Transcripts (3 UTRs)",
                                                            "Abnormal Transcripts (4 UTRs)",
                                                            "Abnormal Transcripts (5 UTRs)",
                                                            "Abnormal Transcripts (6+ UTRs)"))

  return(list(all_ratios,stringtie_UTR2exon_normalUTR_ratio,stringtie_UTR2exon_abnormalUTR_3c_ratio))
}



# #### make plots function - this was not used in final outputs
# 
# plotter <- function(out_list,isWB,title_idx) {
#   ratios <- as.data.frame(out_list[1])
#   normals <- as.data.frame(out_list[2])
#   abnormals_3UTR <- as.data.frame(out_list[3])
#   
#   ratios_counts <- ratios %>% dplyr::group_by(source) %>% 
#     dplyr::mutate(count_lab=paste0("Transcripts: ",as.character(n()))) %>% 
#     dplyr::mutate(count_num=n()) %>% 
#     dplyr::distinct(count_num,.keep_all = T) %>% 
#     dplyr::ungroup()
#   ratios_wP <- ratios_counts %>%
#     dplyr::mutate(percent=round(count_num/ratios_counts$count_num[1]*100,digits = 2)) %>%
#     dplyr::mutate(percent_lab=paste0(as.character(percent),"%"))
#   
#   tr_labels <- data.frame(source2=c("All Transcripts",
#                          "Normal Transcripts (1-2 UTRs)",
#                          "Abnormal Transcripts (3 UTRs)",
#                          "Abnormal Transcripts (4 UTRs)",
#                          "Abnormal Transcripts (5 UTRs)",
#                          "Abnormal Transcripts (6+ UTRs)"), label=ratios_wP$count_lab)
#   
#   pc_labels <- data.frame(source2=c("All Transcripts",
#                                     "Normal Transcripts (1-2 UTRs)",
#                                     "Abnormal Transcripts (3 UTRs)",
#                                     "Abnormal Transcripts (4 UTRs)",
#                                     "Abnormal Transcripts (5 UTRs)",
#                                     "Abnormal Transcripts (6+ UTRs)"), label=ratios_wP$percent_lab)
#   
#   
#   pc_labels$source2 <- factor(pc_labels$source2, levels = c("All Transcripts",
#                                                               "Normal Transcripts (1-2 UTRs)",
#                                                               "Abnormal Transcripts (3 UTRs)",
#                                                               "Abnormal Transcripts (4 UTRs)",
#                                                               "Abnormal Transcripts (5 UTRs)",
#                                                               "Abnormal Transcripts (6+ UTRs)"))
#   tr_labels$source2 <- factor(tr_labels$source2, levels = c("All Transcripts",
#                                                             "Normal Transcripts (1-2 UTRs)",
#                                                             "Abnormal Transcripts (3 UTRs)",
#                                                             "Abnormal Transcripts (4 UTRs)",
#                                                             "Abnormal Transcripts (5 UTRs)",
#                                                             "Abnormal Transcripts (6+ UTRs)"))
#   
#  
#   
#   if(isWB==1) {
#     #ratio hist
#     p1<- ggplot() + geom_histogram(data=ratios,aes(ratio),bins=100) + facet_wrap(~source2) + theme_bw() +ylab("Transcript count")+xlab("CDS length / Transcript length")
#     p1.1 <- p1 + geom_text(x=0.25,y=y_pos[title_idx],aes(label=label),data = tr_labels) + geom_text(x=0.25,y=y_pos[title_idx]-100,aes(label=label),data = pc_labels) +ggtitle(title_vector[title_idx]) +xlim(0,0.999)
#   
#     
#     #normalized ratio hist
#     p2<-ggplot() + geom_histogram(data=ratios,aes(ratio,y = 0.01*..density..),bins=100) + facet_wrap(~source2) + theme_bw() +ylab("Transcript fraction")+xlab("CDS length / Transcript length")
#     p2.1 <- p2 + geom_text(x=0.25,y=0.04,aes(label=label),data = tr_labels) + geom_text(x=0.25,y=0.035,aes(label=label),data = pc_labels)+ggtitle(title_vector[title_idx])+xlim(0,0.999)
# 
#   } else {
#     #ratio hist
#     p1<- ggplot() + geom_histogram(data=ratios,aes(ratio),bins=100) + facet_wrap(~source2) + theme_bw() +ylab("Transcript count")+xlab("CDS length / Transcript length")
#     p1.1 <- p1 + geom_text(x=0.25,y=y_pos[title_idx],aes(label=label),data = tr_labels) + geom_text(x=0.25,y=y_pos[title_idx]-50,aes(label=label),data = pc_labels) +ggtitle(title_vector[title_idx])
#     
#     
#     #normalized ratio hist
#     p2<-ggplot() + geom_histogram(data=ratios,aes(ratio,y = 0.01*..density..),bins=100) + facet_wrap(~source2) + theme_bw() +ylab("Transcript fraction")+xlab("CDS length / Transcript length")
#     p2.1 <- p2 + geom_text(x=0.25,y=0.04,aes(label=label),data = tr_labels) + geom_text(x=0.25,y=0.035,aes(label=label),data = pc_labels)+ggtitle(title_vector[title_idx])
#     
#   }
#   #exon_count vs ratio scatterplot
#   p3 <- ggplot() + geom_point(data=normals,aes(x=exon_count,y=ratio))+ theme_bw() +xlab("Exon Count")+ylab("CDS length / Exon length")+ggtitle(title_vector[title_idx])
# 
#   # exon_sum vs ratio scatterplot
#   #normals_lab <- dplyr::mutate(normals,fill=ifelse(Parent_ID %in% shortCDS_single_exon$Parent_ID | Parent_ID %in% normals_multiExon_single_CDS$Parent_ID,"removed","kept"))
#   p4 <- ggplot() + geom_point(data=normals,aes(y=exon_sum,x=ratio))+ theme_bw() +ylab("Transcript Length")+xlab("CDS length / Exon length") +ggtitle(title_vector[title_idx])
# 
#   #abnormals 3s
#   p5 <- ggplot() + geom_point(data=abnormals_3UTR,aes(x=exon_count,y=ratio))+ theme_bw() +xlab("Exon Count")+ylab("CDS length / Exon length") +ggtitle(title_vector[title_idx])
# 
#   return(list(p1.1,p2.1,p3,p4,p5))
# }

### cleanup GFFs function
cleanup_abnormal4p <- function(gff,normals,abnormals3c) {
  
  normals_shortCDS_single_exon <- normals %>% dplyr::filter(!(exon_count == 1 & ratio < 0.75))
  
  normals_w3c <- bind_rows(normals_shortCDS_single_exon,abnormals3c)
  
  stringtie_L1_features <- gff %>% dplyr::filter(type=="gene") %>% 
    tidyr::separate(attributes,into=c("Feature_ID","Name"),sep=";")
  
  stringtie_L2_features <- gff %>% dplyr::filter(type=="mRNA")%>% 
    tidyr::separate(attributes,into=c("Feature_ID","Parent_ID","Name","orf-model"),sep=";") %>%
    dplyr::group_by(Parent_ID) %>%
    dplyr::mutate(transcript_count=n()) %>%
    dplyr::ungroup()
  stringtie_L2_features$Parent_ID <- gsub(".*=","ID=",stringtie_L2_features$Parent_ID)
  
  stringtie_L3_features <- gff %>% dplyr::filter(grepl("exon|intron|CDS|start_codon|stop_codon|five_prime_UTR|three_prime_UTR",type))%>% 
    tidyr::separate(attributes,into=c("Feature_ID","Parent_ID"),sep=";")
  stringtie_L3_features$Parent_ID <- gsub(".*=","ID=",stringtie_L3_features$Parent_ID)
  
  L1L2_relationships <- stringtie_L2_features %>% dplyr::select(Feature_ID,Parent_ID)
  colnames(L1L2_relationships) <- c("mRNA_ID","Gene_ID")
  

  L2_toRM <- stringtie_L2_features %>% dplyr::filter(!Feature_ID %in% normals_w3c$Parent_ID) %>% 
    group_by(Parent_ID) %>%
    dplyr::mutate(current_transcripts=n()) %>%
    dplyr::ungroup()
  
  L1_toRM <- L2_toRM %>% dplyr::filter(transcript_count==current_transcripts) %>% dplyr::distinct(Parent_ID,.keep_all = T) %>% dplyr::select(Parent_ID)
  
  #filtered_L1_features_split_attributes <- stringtie_L1_features %>% dplyr::filter(!Feature_ID %in% L1_toRM$Parent_ID) 
  filtered_L1_features  <- stringtie_L1_features %>% dplyr::filter(!Feature_ID %in% L1_toRM$Parent_ID) %>%
    tidyr::unite("attributes",Feature_ID,Name,sep=";")
  
  filtered_L2_features <- stringtie_L2_features %>% dplyr::filter(Feature_ID %in% normals_w3c$Parent_ID) 
  filtered_L2_features$Parent_ID <- gsub(".*=","Parent=",filtered_L2_features$Parent_ID)
  filtered_L2_features <- filtered_L2_features %>%  tidyr::unite("attributes",Feature_ID,Parent_ID,Name,`orf-model`,sep=";",na.rm = T) %>% dplyr::select(-transcript_count)
  
  filtered_L3_features <- stringtie_L3_features %>% dplyr::filter(Parent_ID %in% normals_w3c$Parent_ID) 
  filtered_L3_features$Parent_ID <- gsub(".*=","Parent=",filtered_L3_features$Parent_ID)
  filtered_L3_features <- filtered_L3_features %>% tidyr::unite("attributes",Feature_ID,Parent_ID,sep=";",na.rm =T)
  
  filtered_abnormals_GFF <- bind_rows(filtered_L1_features,filtered_L2_features,filtered_L3_features) %>% dplyr::arrange(seqid,start)
  
  dropped_L1_features <- stringtie_L1_features %>% dplyr::filter(Feature_ID %in% L1_toRM$Parent_ID) %>%
    tidyr::unite("attributes",Feature_ID,Name,sep=";")
  
  dropped_L2_features <- stringtie_L2_features %>% dplyr::filter(!Feature_ID %in% normals_w3c$Parent_ID) 
  dropped_L2_features$Parent_ID <- gsub(".*=","Parent=",dropped_L2_features$Parent_ID)
  dropped_L2_features <- dropped_L2_features %>%  tidyr::unite("attributes",Feature_ID,Parent_ID,Name,`orf-model`,sep=";",na.rm = T) %>% dplyr::select(-transcript_count)
  
  dropped_L3_features <- stringtie_L3_features %>% dplyr::filter(!Parent_ID %in% normals_w3c$Parent_ID) 
  dropped_L3_features$Parent_ID <- gsub(".*=","Parent=",dropped_L3_features$Parent_ID)
  dropped_L3_features <- dropped_L3_features %>% tidyr::unite("attributes",Feature_ID,Parent_ID,sep=";",na.rm =T)
  
  dropped_abnormals_GFF <- bind_rows(dropped_L1_features,dropped_L2_features,dropped_L3_features) %>% dplyr::arrange(seqid,start)
  
  return(list(filtered_abnormals_GFF,dropped_abnormals_GFF))
}


#read gff
file <- args[1]
gff <- read.gff(file)

#get dir and filename
dir <- dirname(file)
fileName <- basename(file)

# pull lists of normal transcripts with <2 non-coding exons (<4+ UTR features)
out_list <- pull_ratios(gff,args[2])

# remove abnormal transcripts
filtered_GFFs <- cleanup_abnormal4p(gff,out_list[[2]],out_list[[3]]) 

#set output file names
out1 <- gsub(".gff",".final_kept.gff",fileName)
out2 <- gsub(".gff",".final_dropped.gff",fileName)
outName1 <- paste0(dir,"/",out1)
outName2 <- paste0(dir,"/",out2)

#write GFF of kept and removed transcripts
write.table(filtered_GFFs[[1]],file=outName1,quote = F,sep = "\t",row.names = F,col.names = F,na = ".")
write.table(filtered_GFFs[[2]],file=outName2,quote = F,sep = "\t",row.names = F,col.names = F,na = ".")

