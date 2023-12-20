args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)


#load QX1410 protein length ratios
braker_QX <- read.table(args[1], header = FALSE, sep = "", dec = ".")
stringtie_QX <- read.table(args[2], header = FALSE, sep = "", dec = ".")
merger_QX <- read.table(args[3], header = FALSE, sep = "", dec = ".")

#read VX34 protein legnth ratios
braker_VX <- read.table(args[4], header = FALSE, sep = "", dec = ".")
stringtie_VX <- read.table(args[5], header = FALSE, sep = "", dec = ".")
merger_VX <- read.table(args[6], header = FALSE, sep = "", dec = ".")

#read AF16 protein length ratios
AF16_WS255 <- read.table(args[7], header = FALSE, sep = "", dec = ".")
AF16_WS280 <- read.table(args[8], header = FALSE, sep = "", dec = ".")


colnames(braker_QX) <- c('BRAKER')
colnames(stringtie_QX) <- c('StringTie')
colnames(merger_QX) <- c('Merger')

colnames(braker_VX) <- c('BRAKER')
colnames(stringtie_VX) <- c('StringTie')
colnames(merger_VX) <- c('Merger')

colnames(AF16_WS255) <- c('AF16.WS255')
colnames(AF16_WS280) <- c('AF16.WS280')


braker_QX <- data.frame(as.table(as.matrix(braker_QX)))
braker_QX <- braker_QX[,2:3]
stringtie_QX <- data.frame(as.table(as.matrix(stringtie_QX)))
stringtie_QX <- stringtie_QX[,2:3]
merger_QX <- data.frame(as.table(as.matrix(merger_QX)))
merger_QX <- merger_QX[,2:3]

braker_VX <- data.frame(as.table(as.matrix(braker_VX)))
braker_VX <- braker_VX[,2:3]
stringtie_VX <- data.frame(as.table(as.matrix(stringtie_VX)))
stringtie_VX <- stringtie_VX[,2:3]
merger_VX <- data.frame(as.table(as.matrix(merger_VX)))
merger_VX <- merger_VX[,2:3]


AF16_WS255 <- data.frame(as.table(as.matrix(AF16_WS255)))
AF16_WS255 <- AF16_WS255[,2:3]
AF16_WS280 <- data.frame(as.table(as.matrix(AF16_WS280)))
AF16_WS280 <- AF16_WS280[,2:3]


hb_QX <- hist(braker_QX$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
db_QX <- data.frame(x = hb_QX$breaks,y = c(hb_QX$counts,NA))

hs_QX <- hist(stringtie_QX$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
ds_QX <- data.frame(x = hs_QX$breaks,y = c(hs_QX$counts,NA))

hm_QX <- hist(merger_QX$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
dm_QX <- data.frame(x = hm_QX$breaks,y = c(hm_QX$counts,NA))

hb_VX <- hist(braker_VX$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
db_VX <- data.frame(x = hb_VX$breaks,y = c(hb_VX$counts,NA))

hs_VX <- hist(stringtie_VX$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
ds_VX <- data.frame(x = hs_VX$breaks,y = c(hs_VX$counts,NA))

hm_VX <- hist(merger_VX$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
dm_VX <- data.frame(x = hm_VX$breaks,y = c(hm_VX$counts,NA))


h255_AF16 <- hist(AF16_WS255$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
d255_AF16 <- data.frame(x = h255_AF16$breaks,y = c(h255_AF16$counts,NA))

h280_AF16 <- hist(AF16_WS280$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
d280_AF16 <- data.frame(x = h280_AF16$breaks,y = c(h280_AF16$counts,NA))

dbind <- bind_rows(list(BRAKER=db_QX,StringTie=ds_QX,Merger=dm_QX,AF16.WS255=d255_AF16,AF16.WS280=d280_AF16), .id="id")

dbind$id <- factor(dbind$id, levels=c("AF16.WS255","AF16.WS280","BRAKER","Merger","StringTie")) 

library("ggsci")
bot <- ggplot() + 
  geom_step(data = dbind,aes(x = x,y = y,color=id),stat = "identity") +
  xlim(0,2) + coord_cartesian(ylim = c(0, 1500), xlim = c(0.5, 1.5)) + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text = element_text(size=14),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + scale_colour_manual(values=c("grey","blue","chartreuse2","darkgoldenrod1","red"))

top <- ggplot() + 
  geom_step(data = dbind,aes(x = x,y = y,color=id),stat = "identity") +
  xlim(0,2) + coord_cartesian(ylim = c(1550, 7000), xlim = c(0.5, 1.5)) +
  theme(axis.line.x=element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"), 
        legend.box.just = "right", legend.box.background=element_rect(colour ="white", fill ="grey"),
        legend.margin = margin(6, 6, 6, 6)) + scale_colour_manual(values=c("grey","blue","chartreuse2","darkgoldenrod1","red"))

figureQX <- plot_grid(top,bot,rel_heights = c(1,3),align = 'v',ncol=1)



###VX34

dbind <- bind_rows(list(BRAKER=db_VX,StringTie=ds_VX,Merger=dm_VX,AF16.WS255=d255_AF16,AF16.WS280=d280_AF16), .id="id")

dbind$id <- factor(dbind$id, levels=c("AF16.WS255","AF16.WS280","BRAKER","Merger","StringTie"))

library("ggsci")
bot <- ggplot() +
  geom_step(data = dbind,aes(x = x,y = y,color=id),stat = "identity") +
  xlim(0,2) + coord_cartesian(ylim = c(0, 1400), xlim = c(0.5, 1.5)) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text = element_text(size=14),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + scale_colour_manual(values=c("grey","blue","chartreuse2","darkgoldenrod1","red"))

top <- ggplot() +
  geom_step(data = dbind,aes(x = x,y = y,color=id),stat = "identity") +
  xlim(0,2) + coord_cartesian(ylim = c(1550, 7000), xlim = c(0.5, 1.5)) +
  theme(axis.line.x=element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right", legend.box.background=element_rect(colour ="white", fill ="grey"),
        legend.margin = margin(6, 6, 6, 6)) + scale_colour_manual(values=c("grey","blue","chartreuse2","darkgoldenrod1","red"))

figureVX <- plot_grid(top,bot,rel_heights = c(1,3),align = 'v',ncol=1)

y.grob <- grid::textGrob("counts",
                         gp=gpar( col="black", fontsize=17), rot=90)
x.grob <- grid::textGrob("Predicted protein length / N2 protein length",
                         gp=gpar(col="black", fontsize=17))
lab_A <- grid::textGrob("A", x = 0, hjust = 0, gp=gpar(col="black", fontsize=17,fontface="bold"))
lab_B <- grid::textGrob("B", x = 0, hjust = 0, gp=gpar(col="black", fontsize=17,fontface="bold"))

figure_final <- grid.arrange(arrangeGrob(arrangeGrob(figureQX, top=lab_A),arrangeGrob(figureVX, top=lab_B),ncol=2, left = y.grob, bottom = x.grob))

ggsave("FigS7-def.size.pdf", plot=figure_final,device="pdf")

