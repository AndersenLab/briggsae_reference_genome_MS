args <- commandArgs(trailingOnly = TRUE)

library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(zoo)

file <- args[1]
mergerGFF <- read.gff(file)
mergerGFF_genes <- mergerGFF %>% dplyr::filter(type=="gene") %>% 
  tidyr::separate(attributes,into=c("Gene_ID","Name"),sep=";")

###search for overlaps - iterate until no overlaps are found

### minus stranded overlap search
minusGenes <- mergerGFF_genes %>% dplyr::filter(strand=="-") %>% 
                                  dplyr::arrange(seqid,start,-end) %>% 
                                  dplyr::mutate(overlap=ifelse(lag(seqid)==seqid & lag(end) >= start & lag(end) >= end,"overlap_donor",ifelse(lead(seqid)==seqid & lead(start) <= end & lead(end) <= end,"overlap_receiver","N")))

minusDonors <- minusGenes %>% dplyr::filter(overlap=="overlap_donor") 
minusReceivers <- minusGenes %>% dplyr::filter(overlap=="overlap_receiver") 
remain <- minusGenes %>% dplyr::filter(!overlap=="overlap_donor") %>% dplyr::select(-overlap) 

overlaps_present =1 
while(overlaps_present != 0) {
  remain <- remain %>% dplyr::arrange(seqid,start,-end) %>% 
                       dplyr::mutate(overlap=ifelse(lag(seqid)==seqid & lag(end) >= start & lag(end) >= end,"overlap_donor",ifelse(lead(seqid)==seqid & lead(start) <= end & lead(end) <= end,"overlap_receiver","N")))

  if(nrow(dplyr::filter(remain,overlap=="overlap_donor")) == 0) {
    overlaps_present=0
    next
  } else {
    minusDonors <- bind_rows(minusDonors,dplyr::filter(remain,overlap=="overlap_donor"))
    remain <- remain %>% dplyr::filter(!overlap=="overlap_donor") %>% dplyr::select(-overlap) %>% dplyr::filter(!is.na(seqid))
  }
}

### store fully overlapped minus strand genes
minusGenes_wOverlap <- bind_rows(remain,minusDonors) %>% dplyr::arrange(seqid,start,-end)

### plus stranded overlap search
plusGenes <- mergerGFF_genes %>% dplyr::filter(strand=="+") %>% 
  dplyr::arrange(seqid,start,-end) %>% 
  dplyr::mutate(overlap=ifelse(lag(seqid)==seqid & lag(end) >= start & lag(end) >= end,"overlap_donor",ifelse(lead(seqid)==seqid & lead(start) <= end & lead(end) <= end,"overlap_receiver","N")))

plusDonors <- plusGenes %>% dplyr::filter(overlap=="overlap_donor") 
plusReceivers <- plusGenes %>% dplyr::filter(overlap=="overlap_receiver") 

remain <- plusGenes %>% dplyr::filter(!overlap=="overlap_donor") %>% dplyr::select(-overlap) 

overlaps_present =1 
while(overlaps_present != 0) {
  remain <- remain %>% dplyr::arrange(seqid,start,-end) %>% 
    dplyr::mutate(overlap=ifelse(lag(seqid)==seqid & lag(end) >= start & lag(end) >= end,"overlap_donor",ifelse(lead(seqid)==seqid & lead(start) <= end & lead(end) <= end,"overlap_receiver","N")))
  
  if(nrow(dplyr::filter(remain,overlap=="overlap_donor")) == 0) {
    overlaps_present=0
    next
  } else {
    plusDonors <- bind_rows(plusDonors,dplyr::filter(remain,overlap=="overlap_donor"))
    remain <- remain %>% dplyr::filter(!overlap=="overlap_donor") %>% dplyr::select(-overlap) %>% dplyr::filter(!is.na(seqid))
  }
}

### store fully overlapped plus stranded genes
plusGenes_wOverlap <- bind_rows(remain,plusDonors) %>% dplyr::arrange(seqid,start,-end)

### merge plus and minus
mergerGFF_wOverlap <- bind_rows(plusGenes_wOverlap, minusGenes_wOverlap) %>% dplyr::arrange(seqid,strand,start,-end)

###write list of genes fully overlapped by other genes
all_donors <- bind_rows(plusDonors, minusDonors) %>% dplyr::arrange(seqid,start,-end)


###identify parent-child relationships between genes that overlap (Child will be fused as a transcript of Parent)
newParents <- bind_rows(minusReceivers,plusReceivers) %>% dplyr::arrange(seqid,strand,start,-end,Gene_ID) %>% dplyr::mutate(relationship="Parent")
newParents$group <- 1:nrow(newParents)
newParents$gene_feature <- newParents$Gene_ID
newChildren <- bind_rows(minusDonors,plusDonors) %>% dplyr::arrange(seqid,strand,start,-end,Gene_ID) %>% dplyr::mutate(relationship="Child")

### use na.locf to 'fill' new parents of children that will be fused
Overlap_relationships <- bind_rows(newParents,newChildren)%>% dplyr::arrange(seqid,strand,start,-end,Gene_ID)
Overlap_relationships$group <- na.locf(Overlap_relationships$group) 
Overlap_relationships$gene_feature <- na.locf(Overlap_relationships$gene_feature)
Overlap_relationships <- Overlap_relationships %>% dplyr::select(seqid, source, start, end, strand, Gene_ID, relationship,group,gene_feature)

#get dir and filename
dir <- dirname(file)
fileName <- basename(file)

#set output file names
out1 <- gsub(".gff",".donors.txt",fileName)
out2 <- gsub(".gff",".parents.coordinates",fileName)
outName1 <- paste0(dir,"/",out1)
outName2 <- paste0(dir,"/",out2)

### write parent-child relationships
write.table(all_donors %>% dplyr::select(-overlap,-Name),file=outName1,quote = F,sep = "\t")
write.table(newParents %>% dplyr::select(seqid,start,end),file=outName2,quote = F,sep = "\t",row.names = F,col.names = F)
