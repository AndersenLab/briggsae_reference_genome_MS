library(qtl)
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(linkagemapping)
library(zoo)
library(grid)
library(segmented)
library(ggplotify)
library(cowplot)

#general function to set rownames in a dataframe
setRowNames <- function(df, row.names) {
  rownames(df) <- row.names
  df
}

cr.obj <- read.cross(format ="csvr", file="/Users/nic/Downloads/CbRILgeno-full-reduced.csv", estimate.map = F, genotypes = c("Q", "V"),na.strings = "NA",crosstype="riself")

#drop markers flagged for segregation distortion
gt <- geno.table(cr.obj)
gt[gt$P.value < 0.05/totmar(cr.obj),]
todrop <- rownames(gt[gt$P.value < 1e-10,])
cr.obj <- drop.markers(cr.obj, todrop)

#Study genotype frequencies
# g <- pull.geno(dropped_obj)
# gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
# gfreq <- t(t(gfreq) / colSums(gfreq))
# par(mfrow=c(1,2), las=1)
# for(i in 1:2) {
#    plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "BB")[i],
#           ylim=c(0,1))
# }

#estimate recombination fraction
cr.obj_rf <- est.rf(cr.obj)
#plotRF(cr.obj_rf,alternate.chrid = TRUE)
#rf <- pull.rf(cr.obj_rf)
#lod <- pull.rf(cr.obj_rf, what="lod")
#plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#form linkage groups
lg2 <- formLinkageGroups(cr.obj_rf, max.rf=0.35, min.lod=6,reorgMarkers = FALSE)
table(lg2[,2])

cr.obj_lg <- formLinkageGroups(cr.obj_rf, max.rf=0.35, min.lod=6,reorgMarkers = TRUE)

#remove markers outside six major linkage groups
mn7 <- markernames(cr.obj_lg, chr=7)
cr.obj_rf <- drop.markers(cr.obj_rf,mn7)

#ensure only 6 linkage groups remain
lg <-formLinkageGroups(cr.obj_rf, max.rf=0.35, min.lod=6,reorgMarkers = FALSE)
table(lg[,2])

#plot recombination fraction vs LOD
# par(mfrow=c(2,1)) 
# plot(rf, mn7[1], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
# abline(h=0.5, lty=2)
# plot(lod, mn7[1], bandcol="gray70", alternate.chrid=TRUE)



#select genotype error rate for map estimation
# loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
# for(i in seq(along=err)) {
#      cat(i, "of", length(err), "\n")
#      tempmap <- est.map(cr.obj_rf, error.prob=err[i])
#    loglik[i] <- sum(sapply(tempmap, attr, "loglik")) 
# }
# lod <- (loglik - max(loglik))/log(10)
# plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02), ylab=expression(paste(log[10], " likelihood")))


#estimate genetic map
CBmap <- est.map(cr.obj_rf,error.prob=0.005,verbose=TRUE)
CB_newmap <- replacemap(cr.obj_rf,CBmap)
CBmap2 <- est.map(CB_newmap, err=0.005, maxit=100000,verbose=TRUE)
CB_newmap <- replacemap(cr.obj_rf,CBmap2)


####Detect allele Frequency local deviations
###Generate a set of pairwise changes in allele frequency
chromFreq <- list()
a=1
for(c in CB_newmap$geno) {
  gntyp <- data.frame(c$data)
  freqList <- c()
  print(nrow(gntyp))
  print(ncol(gntyp))
  for(i in 1:ncol(gntyp)) {
    marker <- gntyp[,i]
    count <- table(marker)
    count <- as.vector(unname(count[names(count)==1]))
    freq <- count/nrow(gntyp)
    freqList <- c(freqList,freq)
  }
  freq_df <- data.frame(frequency=as.numeric(freqList), marker=colnames(gntyp)) %>% 
                        tidyr::separate(marker,into=c("chrom","str_pos"),remove=TRUE,sep="_")
  freq_df$str_pos <- as.numeric(freq_df$str_pos)
  chromFreq[[a]] <- freq_df
  a <- a+1
}

### Forward search for local deviations in AF
cutoff=0.04
chrom <- c("I","II","III","IV","V","X")
switchOrder <- c("N","N","N","N","Y","N")
plotList <- list()
a=1
colors <- c("red","blue")
forwardList <- list()
for(chr in chromFreq) {
  
  freqDiff <- chr
  
  if(switchOrder[a]=="Y") {
    i = 1
    fills <- c()
    trail=0
    while(trail==0) {
      if (abs(freqDiff[i,1]-freqDiff[i+1,1]) > cutoff) {
        i <- i+1
        fills <- c(fills,colors[2])
        trail=1
      } else {
        i <- i+1
        fills <- c(fills,colors[2])
      }
    }
    
    while(length(fills)!=(nrow(freqDiff)-1)) {
      if((freqDiff[i,1]-freqDiff[i+1,1] < cutoff) & (freqDiff[i,1]-freqDiff[i+1,1] > (-cutoff)))  {
        switch=0
        i <- i+1
        fills <- c(fills,colors[switch+1])
      } else if (abs(freqDiff[i,1]-freqDiff[i+1,1]) > cutoff) {
        fills <- c(fills,colors[switch+1])
        switch = 1
        hold <- i
        i <- i+1
        consec = 1
        while( switch==1) {
          if ((freqDiff[hold,1]-freqDiff[i,1] < cutoff) & (freqDiff[hold,1]-freqDiff[i,1] > (-cutoff)))  {
            switch=0
            i<-i+1
            fills <- c(fills,colors[switch+1])
          } else if(length(fills)!=(nrow(freqDiff)-1)){
            if(consec==8){
              fills <- c(fills,"badcall")
              i <- i+1
              switch=0
            } else {
              fills <- c(fills,colors[switch+1])
              i <- i+1
              consec <- consec+1
            }
          } else {
            break
          }
        }
      }
    }
    fills <- c(fills, tail(fills,n=1))
    freqDiff$fill2 <- fills
    forwardList[[a]] <- fills
    
    p <- ggplot() +
      geom_point(data=freqDiff,aes(str_pos,frequency,colour=fill2,alpha=0.6)) +
      ggtitle(chrom[a]) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none") +
      xlab("Physical Position (Mb)") + ylim(0,1)+
      ylab("Frequency") #+ scale_x_continuous(labels = unit_format(unit="",scale = 1e-6))
    plotList[[a]] <- p
    a <- a+1
  } else {
    i = 1
    fills <- c()
    while(length(fills)!=(nrow(freqDiff)-1)) {
      
      if((freqDiff[i,1]-freqDiff[i+1,1] < cutoff) & (freqDiff[i,1]-freqDiff[i+1,1] > (-cutoff)))  {
        switch=0
        i <- i+1
        fills <- c(fills,colors[switch+1])
      } else if (abs(freqDiff[i,1]-freqDiff[i+1,1]) > cutoff) {
        fills <- c(fills,colors[switch+1])
        switch = 1
        hold <- i
        i <- i+1
        consec=1
        while( switch==1) {
          if ((freqDiff[hold,1]-freqDiff[i,1] < cutoff) & (freqDiff[hold,1]-freqDiff[i,1] > (-cutoff)))  {
            switch=0
            i<-i+1
            fills <- c(fills,colors[switch+1])
          } else if(length(fills)!=(nrow(freqDiff)-1)){
            if(consec==8 & !(a==1)){
              fills <- c(fills,"badcall")
              i <- i+1
              switch=0
            } else {
              fills <- c(fills,colors[switch+1])
              i <- i+1
              consec <- consec+1
            }
          } else {
            break
          }
        }
      }
    }
    fills <- c(fills, tail(fills,n=1))
    freqDiff$fill2 <- fills
    forwardList[[a]] <- fills
    
    p <- ggplot() +
      geom_point(data=freqDiff,aes(str_pos,frequency,colour=fill2,alpha=0.6)) +
      ggtitle(chrom[a]) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none") +
      xlab("Physical Position (Mb)") + ylim(0,1) +
      ylab("Frequency") #+ scale_x_continuous(labels = unit_format(unit="",scale = 1e-6))
    plotList[[a]] <- p
    a <- a+1
  }
  
}
p <- grid.arrange(grobs=plotList,ncol=3)

### Backward search for local deviations in AF
a=1
colors <- c("red","blue")
backList <- list()
for(chr in chromFreq) {
  
  freqDiff <- chr
  
  if(switchOrder[a]=="Y") {
    i = nrow(freqDiff)
    fills <- c()
    
    while(length(fills)!=(nrow(freqDiff)-1)) {
      
      if((freqDiff[i,1]-freqDiff[i-1,1] < cutoff) & (freqDiff[i,1]-freqDiff[i-1,1] > (-cutoff)))  {
        switch=0
        i <- i-1
        fills <- c(fills,colors[switch+1])
      } else if (abs(freqDiff[i,1]-freqDiff[i-1,1]) > cutoff) {
        fills <- c(fills,colors[switch+1])
        switch = 1
        hold <- i
        i <- i-1
        consec=1
        while( switch==1) {
          if ((freqDiff[hold,1]-freqDiff[i,1] < cutoff) & (freqDiff[hold,1]-freqDiff[i,1] > (-cutoff)))  {
            switch=0
            i<-i-1
            fills <- c(fills,colors[switch+1 ])
          } else if(length(fills)!=(nrow(freqDiff)-1)){
            if(consec==8){
              fills <- c(fills,"badcall")
              i <- i-1
              switch=0
            } else {
              fills <- c(fills,colors[switch+1])
              i <- i-1
              consec <- consec+1
            }
          } else {
            break
          }
        }
      }
    }
    fills <- c(fills, tail(fills,n=1))
    freqDiff$fill2 <- rev(fills)
    backList[[a]] <- rev(fills)
    
    p <- ggplot() +
      geom_point(data=freqDiff,aes(str_pos,frequency,colour=fill2,alpha=0.6)) +
      ggtitle(chrom[a]) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none") +
      xlab("Physical Position (Mb)") + ylim(0,1)+
      ylab("Frequency") #+ scale_x_continuous(labels = unit_format(unit="",scale = 1e-6))
    plotList[[a+6]] <- p
    a <- a+1
  } else {
    i = nrow(freqDiff)
    fills <- c()
    while(length(fills)!=(nrow(freqDiff)-1)) {
      
      if((freqDiff[i,1]-freqDiff[i-1,1] < cutoff) & (freqDiff[i,1]-freqDiff[i-1,1] > (-cutoff)))  {
        switch=0
        i <- i-1
        fills <- c(fills,colors[switch+1])
      } else if (abs(freqDiff[i,1]-freqDiff[i-1,1]) > cutoff) {
        fills <- c(fills,colors[switch+1])
        switch = 1
        hold <- i
        i <- i-1
        consec = 1
        while( switch==1) {
          
          if ((freqDiff[hold,1]-freqDiff[i,1] < cutoff) & (freqDiff[hold,1]-freqDiff[i,1] > (-cutoff)))  {
            switch=0
            i<-i-1
            fills <- c(fills,colors[switch+1])
          } else if(length(fills)!=(nrow(freqDiff)-1)){
            if(consec==8){
              fills <- c(fills,"badcall")
              i <- i-1
              switch=0
            } else {
              fills <- c(fills,colors[switch+1])
              i <- i-1
              consec <- consec+1
            }
          } else {
            break
          }
        }
      }
    }
    fills <- c(fills, tail(fills,n=1))
    freqDiff$fill2 <- rev(fills)
    backList[[a]] <- rev(fills)
    
    
    p <- ggplot() +
      geom_point(data=freqDiff,aes(str_pos,frequency,colour=fill2,alpha=0.6)) +
      ggtitle(chrom[a]) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none") +
      xlab("Physical Position (Mb)") + ylim(0,1) +
      ylab("Frequency") #+ scale_x_continuous(labels = unit_format(unit="",scale = 1e-6))
    plotList[[a+6]] <- p
    a <- a+1
  }
  
}
p <- grid.arrange(grobs=plotList,ncol=3)


a=1
allflags <- list()
clean_set <- list()
for(chr in chromFreq) {
  
  flags <- data.frame(forwardList[a],backList[a])
  colnames(flags) <- c("forward","backward")
  frqFlags <- cbind(chr,flags)
  frqFlags <- frqFlags %>% mutate(flag=ifelse(forward=="red" & backward=="red","red",ifelse(forward=="blue" & backward=="blue","blue","red")))
  
  if(a==2) {
    p <- ggplot(frqFlags,aes(str_pos,frequency,colour="blue",alpha=0.6)) +
      geom_point() +
      ggtitle(chrom[a]) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none") + scale_color_manual(values="blue") +
      xlab("Physical Position (Mb)") + ylim(0,1) +
      ylab("Frequency") #+ scale_x_continuous(labels = unit_format(unit="",scale = 1e-6))
    plotList[[a+12]] <- p
    allflags[[a]] <- frqFlags
    clean_set[[a]] <- frqFlags[frqFlags$flag=="red",]
    a<-a+1
  } else {
    p <- ggplot() +
      geom_point(data=frqFlags,aes(str_pos,frequency,colour=flag,alpha=0.6)) +
      ggtitle(chrom[a]) + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none") + scale_color_manual(values=c("red","blue")) +
      xlab("Physical Position (Mb)") + ylim(0,1) +
      ylab("Frequency") #+ scale_x_continuous(labels = unit_format(unit="",scale = 1e-6))
    plotList[[a+12]] <- p
    allflags[[a]] <- frqFlags
    clean_set[[a]] <- frqFlags[frqFlags$flag=="red",]
    a<-a+1
  }
}
p <- grid.arrange(grobs=plotList,ncol=3)


all_clean_flags <- list()
a=1
for(c in allflags) {
  clean_fflags <- c()
  clean_bflags <- c()
  fflags <- rev(c$forward)
  bflags <- c$backward
  idx=1
  while(!(idx-1 ==length(bflags))) {
    if (bflags[idx] == "badcall") {
      clean_bflags <- c(clean_bflags,"badcall")
      idx <- idx+1
      while(bflags[idx] == "blue") {
        clean_bflags <- c(clean_bflags,"badcall")
        idx <- idx+1
      }
    } else {
      clean_bflags <- c(clean_bflags,as.character(bflags[idx]))
      idx <- idx+1
    }
  }
  idx=1
  while(!(idx-1 ==length(fflags))) {
    if (fflags[idx] == "badcall") {
      clean_fflags <- c(clean_fflags,"badcall")
      idx <- idx+1
      while(fflags[idx] == "blue") {
        clean_fflags <- c(clean_fflags,"badcall")
        idx <- idx+1
      }
    } else {
      clean_fflags <- c(clean_fflags,as.character(fflags[idx]))
      idx <- idx+1
    }
  }
  chromFlags <- data.frame(rev(clean_fflags),clean_bflags)
  all_clean_flags[[a]] <- chromFlags
  a<-a+1
}


a=1
b=1
allflags <- list()
plotList <- list() 
clean_set <- list()
for(chr in chromFreq) {
  
  flags <- data.frame(forwardList[a],backList[a])
  colnames(flags) <- c("forward","backward")
  frqFlags <- cbind(chr,flags)
  cleanflags <- all_clean_flags[[a]]
  colnames(cleanflags) <- c("clean_fflags","clean_bflags")
  frqFlags <- cbind(frqFlags, cleanflags)
  frqFlags <- frqFlags %>% mutate(flag=ifelse(clean_fflags=="blue" | clean_bflags=="blue","blue","red"))
  
  p <- ggplot() +
      geom_point(data=frqFlags,aes(str_pos,frequency,colour=clean_fflags,alpha=0.6)) +
      ggtitle(chrom[a]) + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none") +# scale_color_manual(values=c("red","blue")) +
      xlab("Physical Position (Mb)") + ylim(0,1) +
      ylab("Frequency") #+ scale_x_continuous(labels = unit_format(unit="",scale = 1e-6))
    
  q <- ggplot() +
    geom_point(data=frqFlags,aes(str_pos,frequency,colour=clean_bflags,alpha=0.6)) +
    ggtitle(chrom[a]) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") +# scale_color_manual(values=c("red","blue")) +
    xlab("Physical Position (Mb)") + ylim(0,1) +
    ylab("Frequency")

  r <- ggplot() +
    geom_point(data=frqFlags,aes(str_pos,frequency,colour=flag,alpha=0.6)) +
    ggtitle(chrom[a]) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") +# scale_color_manual(values=c("red","blue")) +
    xlab("Physical Position (Mb)") + ylim(0,1) +
    ylab("Frequency")

  
    plotList[[b]] <- p
    plotList[[b+1]] <- q
    plotList[[b+2]] <- r
    allflags[[a]] <- frqFlags
    clean_set[[a]] <- frqFlags[frqFlags$flag=="red",]
    a<-a+1
    b<-b+3
}
p <- grid.arrange(grobs=plotList,ncol=3)


#plotList[[15]] + xlim(3000000,5000000) 







flaggedMarkers <- list()
a=1
for(chrFrq in allflags) {
  flaggedMarkers[[a]] <- chrFrq[chrFrq$flag=="blue",]
  a<-a+1
}

#map <- pull.map()

plotList2 <- list()
a=1
for(c in CBmap2) {
  df <- data.frame(names(c),as.vector(c))
  colnames(df) <- c("name","est_map")
  df2 <- df %>% 
    tidyr::separate(name,into=c("chrom","str_pos"),remove=FALSE,sep="_") %>% 
    dplyr::mutate(gen_pos=(as.numeric(est_map-est_map[1]))) %>%
    dplyr::mutate(phys_pos=(as.numeric(str_pos)))
  if(nrow(flaggedMarkers[[a]])==0) {
    df2 <- df2 %>% dplyr::mutate(flag="Y")
  } else {
    df2 <- df2 %>% dplyr::mutate(flag=ifelse(str_pos %in% flaggedMarkers[[a]]$str_pos,"Y","N"))
  }

  p <- ggplot() + 
    geom_point(data=df2,aes(phys_pos,gen_pos)) + 
    ggtitle(chrom[a]) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") + 
    xlab("Physical Position (Mb)") + 
    ylab("Genetic Position (cM)") + scale_x_continuous(labels = unit_format(unit="",scale = 1e-6))
  plotList2[[a]] <- p
  a <- a+1
}
p <- grid.arrange(grobs=plotList,ncol=3)


### Drop markers flagged for local AF deviations
todrop <- bind_rows(flaggedMarkers, .id = "chrID") %>% tidyr::unite(dAF,chrom,str_pos,sep="_",remove=FALSE)
dAF <- todrop$dAF

### Re-estimate genetic map
CB_newmap_dAF <- drop.markers(CB_newmap,dAF)
dAFmap <- est.map(CB_newmap_dAF,error.prob=0.005,verbose=TRUE)
CB_newmap_dAF <- replacemap(CB_newmap_dAF,dAFmap)

library(ggrepel)
plotList2 <- list()
a=1
for(c in dAFmap) {
  df <- data.frame(names(c),as.vector(c))
  colnames(df) <- c("name","est_map")
  df2 <- df %>% 
    tidyr::separate(name,into=c("chrom","str_pos"),remove=FALSE,sep="_") %>% 
    dplyr::mutate(gen_pos=(as.numeric(est_map-est_map[1]))) %>%
    dplyr::mutate(phys_pos=(as.numeric(str_pos)))
  p <- ggplot(data=df2,aes(phys_pos,gen_pos)) + 
    geom_point() + 
    ggtitle(chrom[a]) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") + 
    xlab("Physical Position (Mb)") + 
    ylab("Genetic Position (cM)") + scale_x_continuous(labels = unit_format(unit="",scale = 1e-6)) +
    geom_label_repel(aes(label=phys_pos),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50')
  plotList2[[a]] <- p
  a <- a+1
}
p <- grid.arrange(grobs=plotList2,ncol=3)


plotList3 <- list()
a=1
for(c in dAFmap) {
  df <- data.frame(names(c),as.vector(c))
  colnames(df) <- c("name","est_map")
  df2 <- df %>% 
    tidyr::separate(name,into=c("chrom","str_pos"),remove=FALSE,sep="_") %>% 
    dplyr::mutate(gen_pos=(as.numeric(est_map-est_map[1]))) %>%
    dplyr::mutate(phys_pos=(as.numeric(str_pos)))
  
  p <- ggplot(data=df2,aes(phys_pos,gen_pos)) + 
    geom_point() + 
    ggtitle(chrom[a]) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") + 
    xlab("Physical Position (Mb)") + 
    ylab("Genetic Position (cM)") + scale_x_continuous(labels = unit_format(unit="",scale = 1e-6)) #+
  plotList3[[a]] <- p
  a <- a+1
}
p <- grid.arrange(grobs=plotList3,ncol=3)


### Remove manually curated markers
manually_curated <- c("II_12285","II_2216118","II_2218699","II_16571881","II_16571481","II_11779","II_11398","II_11429","III_9518491","III_9517915","III_9518404","III_9516439","IV_49449","IV_50840","IV_46473","IV_49267","V_14255854","V_14255766","V_14255878","V_14256052","X_21407010","X_21406898")

### Re-estiamte map
CB_newmap_manual <- drop.markers(CB_newmap_dAF,manually_curated)
manual_map <- est.map(CB_newmap_manual,error.prob=0.005,verbose=TRUE)
CB_newmap_manual <- replacemap(CB_newmap_manual,manual_map)


### Plot Marey Maps with genome gaps and segmented linear fits
plotList4 <- list()
breakpoints <- list()
slopes <- list()
ks=2
a=1
for(c in manual_map) {
  df <- data.frame(names(c),as.vector(c))
  colnames(df) <- c("name","est_map")
  df2 <- df %>% 
    tidyr::separate(name,into=c("chrom","str_pos"),remove=FALSE,sep="_") %>% 
    dplyr::mutate(gen_pos=(as.numeric(est_map-est_map[1]))) %>%
    dplyr::mutate(phys_pos=(as.numeric(str_pos)))
  
  
  if (a==2) {
    ####
    dfFP <- df2[df2$phys_pos < 15000000,]
    my.lmFP <- lm(gen_pos ~ phys_pos, data =dfFP)
    
    my.segFP <- segmented(my.lmFP, seg.Z=~phys_pos,npsi=ks)
    my.fittedFP <- fitted(my.segFP)
    my.modelFP <- data.frame(PhysicalDistance = dfFP$phys_pos, RecombinationRate = my.fittedFP)
    
    
    my.slopesFP <- coef(my.segFP)
    b0 <- my.slopesFP[[1]]
    b1 <- my.slopesFP[[2]]
    b2 <- my.slopesFP[[3]]
    b3 <- my.slopesFP[[4]]
    
    
    c1<-b1+b2
    break1<- my.segFP$psi[[3]]
    
    c0 <- b0+b1*break1-c1*break1
    
    d1<- b3+c1
    break2 <- my.segFP$psi[[4]]
    
    d0 <- c0 + c1 * break2 - d1*break2
    breakpoints[[a]] <- setRowNames(as.data.frame(my.segFP$psi)[,2:3],c("L.Arm - Center","Center - R.Arm"))
    

    
  } else {
    my.lm <- lm(gen_pos ~ phys_pos, data =df2)
    my.coef <- coef(my.lm)
    my.seg <- segmented(my.lm, seg.Z=~phys_pos,npsi=ks)
    breakpoints[[a]] <- setRowNames(as.data.frame(my.seg$psi)[,2:3],c("L.Arm - Center","Center - R.Arm"))
    my.fitted <- fitted(my.seg)
    my.model <- data.frame(PhysicalDistance = df2$phys_pos, RecombinationRate = my.fitted)
    my.slopes <- coef(my.seg)
    b0 <- my.slopes[[1]]
    b1 <- my.slopes[[2]]
    b2 <- my.slopes[[3]]
    b3 <- my.slopes[[4]]
    
    
    c1<-b1+b2
    break1<- my.seg$psi[[3]]
    
    c0 <- b0+b1*break1-c1*break1
    
    d1<- b3+c1
    break2 <- my.seg$psi[[4]]
    
    d0 <- c0 + c1 * break2 - d1*break2
    
  }
  
  slopes[[a]] <- c(b1*1000000,c1*1000000,d1*1000000)
  
  if(a==2) {
    p <- ggplot(data=df2,aes(phys_pos,gen_pos)) + 
      geom_point() + 
      ggtitle(chrom[a]) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            text=element_text(size=15),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(margin = margin(t = 10, b = -20),hjust=0.02,size=16),
            legend.position = "none") + 
      scale_x_continuous(labels = unit_format(unit="",scale = 1e-6,accuracy=1)) +
      scale_y_continuous(breaks = c(seq(0,round(max(df2$gen_pos),digits=-1),10))) +
      geom_line(data=my.modelFP, aes(x = PhysicalDistance, y = RecombinationRate)) +
      geom_vline(xintercept=14716998,color="#CC79A7",linetype='dashed')+
      geom_abline(intercept = b0, slope = b1, 
                  aes(colour = "Left arm")) +
      geom_abline(intercept = c0, slope = c1, 
                  aes(colour = "Center"))+
      geom_abline(intercept = d0, slope = d1, 
                  aes(colour = "Right arm"))
    
  } else {
    p <- ggplot(data=df2,aes(phys_pos,gen_pos)) + 
      geom_point() + 
      ggtitle(chrom[a]) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            text=element_text(size=15),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(margin = margin(t = 10, b = -20),hjust=0.02,size = 16),
            legend.position = "none") + 
      scale_x_continuous(labels = unit_format(unit="",scale = 1e-6,accuracy=1)) +
      scale_y_continuous(breaks = c(seq(0,round(max(df2$gen_pos),digits=-1),10))) +
      geom_abline(intercept = b0, slope = b1, 
                  aes(colour = "Left arm")) +
      geom_abline(intercept = c0, slope = c1, 
                  aes(colour = "Center"))+
      geom_abline(intercept = d0, slope = d1, 
                  aes(colour = "Right arm"))
    if(a==1) {
      p <- p + geom_vline(xintercept=9640343,color="#CC79A7",linetype='dashed')
    } else if(a==6) {
      gapx <- c(6344845,13120827,16639414,17893975,22154861)
      p <- p+ geom_vline(xintercept=gapx, color="#CC79A7",linetype='dashed')
    }
  }
          
  plotList4[[a]] <- p
  a <- a+1
}


chromFreq <- list()
a=1
for(c in CB_newmap_manual$geno) {
  gntyp <- data.frame(c$data)
  freqList <- c()
  for(i in 1:ncol(gntyp)) {
    marker <- gntyp[,i]
    count <- table(marker)
    count <- as.vector(unname(count[names(count)==1]))
    freq <- count/nrow(gntyp)
    freqList <- c(freqList,freq)
  }
  freq_df <- data.frame(frequency=as.numeric(freqList), marker=colnames(gntyp)) %>% 
    tidyr::separate(marker,into=c("chrom","str_pos"),remove=TRUE,sep="_")
  freq_df$str_pos <- as.numeric(freq_df$str_pos)
  chromFreq[[a]] <- freq_df
  a <- a+1
}

###Define sliding window and step sizes
windowStart <- seq(0,23000000,5000)
windowEnd <- seq(50000,23050000,5000)
windows <- data.frame(windowStart,windowEnd) 


### Average allele frequencies across sliding window
a=1
allMeans <- list()
for(chr in chromFreq) {
  meanList <- c()
  for(row in 1:nrow(windows)) {
    interval <- windows[row,]
    ivl_markers <- chr %>% dplyr::filter(str_pos > interval$windowStart & str_pos < interval$windowEnd)
    meanFreq <- mean(ivl_markers$frequency)
    meanList <- c(meanList,meanFreq)
  }
  allMeans[[a]] <- meanList
  a<-a+1
}

allWindows <- list()
a=1
for (means in allMeans) {
  test <- cbind(windows,means)
  colnames(test) <- c("windowStart","windowEnd","meanFreq")
  test$meanFreq[is.nan(test$meanFreq)] <- NA
  test <- na.locf(test)
  markers <- data.frame(clean_set[a])
  test <- test %>% dplyr::mutate(chr=chrom[a]) %>%
    dplyr::filter(windowStart<max(markers$str_pos))
  allWindows[[a]] <- test
  a<-a+1
}

finalSet <- bind_rows(allWindows,.id="CHROM.ID")
### Plot allele frequencies
af <- ggplot(finalSet, aes(windowStart,meanFreq,color=chr)) + 
  geom_line() +
  ylim(0,1) +
  geom_hline(yintercept =0.5,color="black",alpha=1,linetype='dashed') + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.95,0.76),
        legend.key = element_rect(fill="white"),
        legend.key.size= unit(0.6,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12,hjust=0.35),
        axis.title = element_text(size=16),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  xlab("Physical position (Mb)") +
  ylab("QX1410 allele frequency") +
  scale_x_continuous(labels = unit_format(unit="",scale = 1e-6,accuracy=1))


### Plot main display figure
plot(manual_map)
marey <- grid.arrange(grobs=plotList4,ncol=3,left=textGrob("Genetic position (cM)",rot=90,gp=gpar(fontsize=16)),bottom=textGrob("Physical position (Mb)", gp=gpar(fontsize=16)))
topseq <- c(rep(1,7),rep(2,993))
botseq <- rep(3,1000)
lay=rbind(topseq,botseq)
af2 <- arrangeGrob(af,top=textGrob("B",gp=gpar(fontsize=20),x = 0, hjust = 0))
grid.arrange(nullGrob(),marey,af2,layout_matrix=lay,top=textGrob("A",gp=gpar(fontsize=20),x = 0, hjust = 0))
##ggsave


### Write genetic distances
map_table <- pull.map(CB_newmap_manual,as.table = T)
write.table(map_table,file="c_briggsae_genetic_distances.txt",quote = F)

gt <- geno.table(CB_newmap_manual, scanone.output=TRUE)
par(mfrow=c(1,1))
plot(gt, ylab=expression(paste(-log[10], " P-value")))
plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")


### Plot rugplot of genetic map
plotMap(CB_newmap_manual, show.marker.names=FALSE)

### Generate table with breakpoint calls with standard error
all_BP <- do.call("rbind",breakpoints)
all_BP$Breakpoint <- rownames(all_BP)
all_BP$Chromosome <- c("I","I","II","II","III","III","IV","IV","V","V","X","X")
all_BP <- all_BP[,c(4,3,1,2)]
write.table(all_BP,file="c_briggsae_arm-center_breakpoints.txt",quote = F,row.names = F,sep="\t")

save(CB_newmap_manual, file="CB_genetic_map.Rda")

### Manually defined tips
xltip <- "1443949"
xrtip <- "22160685"
vltip <- "468664"
vrtip <- "19546908"
ivltip <-"702022"
ivrtip <- "16366490"
iiiltip <- "716441"
iiirtip <- "13994534"
iiltip <- "566445"
iirtip <- "16571409"
iltip <- "388131"
irtip <- "14969701"

tips <- c(iltip,irtip,iiltip,iirtip,iiiltip,iiirtip,ivltip,ivrtip,vltip,vrtip,xltip,xrtip)

plotList5 <- list()
i=1
a=1
for (plot in plotList4){
  p <- plot + geom_vline(xintercept =as.numeric(tips[i]),color='red') +geom_vline(xintercept =as.numeric(tips[i+1]),color='red')
  plotList5[[a]] <- p
  a <- a+1
  i <- i+2
}
mareytips <- grid.arrange(grobs=plotList5,ncol=3,left=textGrob("Genetic position (cM)",rot=90,gp=gpar(fontsize=16)),bottom=textGrob("Physical position (Mb)", gp=gpar(fontsize=16)))



### Plot supplementary figured of flagged markers with local AF deviations
plotList6 <- list()
plotList6[[1]] <- plotList[[3]] + theme(axis.title=element_blank(), plot.title = element_text(margin = margin(t = 10, b = -20),hjust=0.02,size = 16),panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_x_continuous(labels = unit_format(unit="",scale = 1e-6,accuracy=1))
plotList6[[2]] <- plotList[[6]]+ theme(axis.title=element_blank(), plot.title = element_text(margin = margin(t = 10, b = -20),hjust=0.02,size = 16), panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_x_continuous(labels = unit_format(unit="",scale = 1e-6,accuracy=1))
plotList6[[3]] <- plotList[[9]]+ theme(axis.title=element_blank(), plot.title = element_text(margin = margin(t = 10, b = -20),hjust=0.02,size = 16),panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_x_continuous(labels = unit_format(unit="",scale = 1e-6,accuracy=1))
plotList6[[4]] <- plotList[[12]]+ theme(axis.title=element_blank(), plot.title = element_text(margin = margin(t = 10, b = -20),hjust=0.02,size = 16),panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_x_continuous(labels = unit_format(unit="",scale = 1e-6,accuracy=1))
plotList6[[5]] <- plotList[[15]]+ theme(axis.title=element_blank(), plot.title = element_text(margin = margin(t = 10, b = -20),hjust=0.02,size = 16),panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_x_continuous(labels = unit_format(unit="",scale = 1e-6,accuracy=1))
plotList6[[6]] <- plotList[[18]]+ theme(axis.title=element_blank(), plot.title = element_text(margin = margin(t = 10, b = -20),hjust=0.02,size = 16),panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_x_continuous(labels = unit_format(unit="",scale = 1e-6,accuracy=1))
AFskews <- grid.arrange(grobs=plotList6,ncol=3,left=textGrob("QX1410 allele frequency",rot=90,gp=gpar(fontsize=16)),bottom=textGrob("Physical position (Mb)", gp=gpar(fontsize=16)))
plotList7 <- list()
i=1
a=1
for (plot in plotList4){
  p <- plot + geom_vline(xintercept =as.numeric(all_BP$Est.[i]),color='red') +geom_vline(xintercept =as.numeric(all_BP$Est.[i+1]),color='red')
  plotList7[[a]] <- p
  a <- a+1
  i <- i+2
}
arm_center <- grid.arrange(grobs=plotList7,ncol=3,left=textGrob("Genetic position (cM)",rot=90,gp=gpar(fontsize=16)),bottom=textGrob("Physical position (Mb)", gp=gpar(fontsize=16)))

write.csv(matrix(dAF, nrow=1), file ="~/Desktop/dAF_markers.csv", row.names=FALSE,quote = FALSE)

