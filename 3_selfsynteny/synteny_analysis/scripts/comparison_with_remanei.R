#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# sp1
sp1.synteny <- read.table("CBRIG_CREMA.syn", col.names = c("chromosome", "syn", "break", "perc"))
sp1.synteny$chromosome <- factor(sp1.synteny$chromosome)
sp1.synteny.df <- gather(sp1.synteny, X, prop, perc, factor_key=TRUE)
sp1.total <- filter(sp1.synteny.df, chromosome == "Total")
sp1.total <- sp1.total$prop
sp1.synteny.df <- filter(sp1.synteny.df, chromosome != "Total")
sp1.synteny.df$reproductive_mode <- (rep("selfer", times=6))

# sp2
sp2.synteny <- read.table("CNIGO_CREMA.syn", col.names = c("chromosome", "syn", "break", "perc"))
sp2.synteny$chromosome <- factor(sp2.synteny$chromosome)
sp2.synteny.df <- gather(sp2.synteny, X, prop, perc, factor_key=TRUE)
sp2.total <- filter(sp2.synteny.df, chromosome == "Total")
sp2.total <- sp2.total$prop
sp2.synteny.df <- filter(sp2.synteny.df, chromosome != "Total")
sp2.synteny.df$reproductive_mode <- (rep("outcrosser", times=6))

# join
df1 <- rbind(sp1.synteny.df, sp2.synteny.df)

p1 <- ggplot(data=df1, aes(x=chromosome, y=prop, color=reproductive_mode)) + 
  geom_line(color="black") + 
  geom_point(size=4) +
  theme_bw() + 
  scale_color_manual(values = c("#ff8800", "#0092c1")) + 
  ggtitle(expression(paste(italic("C. briggsae/C. nigoni"), " vs ", italic("C. remanei")))) + 
  xlab("Chromosome") + 
  ylab(expression(paste("% neighbouring gene pairs colinear in ", italic("C. remanei")))) + 
  labs(color = "Reproductive mode") + 
  theme(plot.title = element_text(face="bold", size=18), axis.text = element_text(size=12), axis.title = element_text(size=14)) 

# sp3
sp3.synteny <- read.table("CELEG_CREMA.syn", col.names = c("chromosome", "syn", "break", "perc"))
sp3.synteny$chromosome <- factor(sp3.synteny$chromosome)
sp3.synteny.df <- gather(sp3.synteny, X, prop, perc, factor_key=TRUE)
sp3.total <- filter(sp3.synteny.df, chromosome == "Total")
sp3.total <- sp3.total$prop
sp3.synteny.df <- filter(sp3.synteny.df, chromosome != "Total")
sp3.synteny.df$reproductive_mode <- (rep("selfer", times=6))

# sp4
sp4.synteny <- read.table("CINOP_CREMA.syn", col.names = c("chromosome", "syn", "break", "perc"))
sp4.synteny$chromosome <- factor(sp4.synteny$chromosome)
sp4.synteny.df <- gather(sp4.synteny, X, prop, perc, factor_key=TRUE)
sp4.total <- filter(sp4.synteny.df, chromosome == "Total")
sp4.total <- sp4.total$prop
sp4.synteny.df <- filter(sp4.synteny.df, chromosome != "Total")
sp4.synteny.df$reproductive_mode <- (rep("outcrosser", times=6))

# join
df2 <- rbind(sp3.synteny.df, sp4.synteny.df)

p2 <- ggplot(data=df2, aes(x=chromosome, y=prop, color=reproductive_mode)) + 
  geom_line(color="black") + 
  geom_point(size=4) +
  theme_bw() + 
  scale_color_manual(values = c("#ff8800", "#0092c1")) + 
  ggtitle(expression(paste(italic("C. elegans/C. inopinata"), " vs ", italic("C. remanei")))) + 
  xlab("Chromosome") + 
  ylab(expression(paste("% neighbouring gene pairs colinear in ", italic("C. remanei")))) + 
  labs(color = "Reproductive mode") + 
  theme(plot.title = element_text(face="bold", size=18), axis.text = element_text(size=12), axis.title = element_text(size=14)) 

p <- p1 / p2 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20, face="bold"))

ggsave("comparison_with_remanei.pdf", p, width=9, height=10)
ggsave("comparison_with_remanei.png", p, width=9, height=10)
