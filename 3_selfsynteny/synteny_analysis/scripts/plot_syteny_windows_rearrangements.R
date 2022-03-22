library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggregplot)

## CELEG and CBRIG synteny windows
syn_windows <- read.table("CELEG_CBRIG.syn_windows", sep="\t", col.names=c("chr", "start", "end", "syn", "break", "prop"))
syn_windows <- syn_windows %>% mutate(mid = (start + end)/2)

I <- syn_windows %>% filter(chr == "I") %>% 
  ggplot(aes(x=mid, y=prop)) + 
  geom_point() + 
  geom_smooth(color="black", fill="light grey", span=0.75) +
  geom_vline(xintercept=3858*1e3, linetype=2) + geom_vline(xintercept=11040*1e3, linetype=2) + 
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 15072434)) + scale_y_continuous(labels=function(x)x*1e2, limits=c(0,1.05)) +
  facet_grid(~chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

II <- syn_windows %>% filter(chr == "II") %>% 
  ggplot(aes(x=mid, y=prop)) + 
  geom_point() + 
  geom_smooth(color="black", fill="light grey", span=0.75) +
  geom_vline(xintercept=4879*1e3, linetype=2) + geom_vline(xintercept=12020*1e3, linetype=2) + 
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 15279421)) + scale_y_continuous(labels=function(x)x*1e2, limits=c(0,1.05)) +
  facet_grid(~chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

III <- syn_windows %>% filter(chr == "III") %>% 
  ggplot(aes(x=mid, y=prop)) + 
  geom_point() + 
  geom_smooth(color="black", fill="light grey", span=0.75) +
  geom_vline(xintercept=3722*1e3, linetype=2) + geom_vline(xintercept=10340*1e3, linetype=2) + 
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 13783801)) + scale_y_continuous(labels=function(x)x*1e2, limits=c(0,1.05)) +
  facet_grid(~chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

IV <- syn_windows %>% filter(chr == "IV") %>%  
  ggplot(aes(x=mid, y=prop)) + 
  geom_point() + 
  geom_smooth(color="black", fill="light grey", span=0.75) +
  geom_vline(xintercept=3896*1e3, linetype=2) + geom_vline(xintercept=12970*1e3, linetype=2) +
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 17493829)) + scale_y_continuous(labels=function(x)x*1e2, limits=c(0,1.05)) +
  facet_grid(~chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

V <- syn_windows %>% filter(chr == "V") %>% 
  ggplot(aes(x=mid, y=prop)) + 
  geom_point() + 
  geom_smooth(color="black", fill="light grey", span=0.75) +
  geom_vline(xintercept=5897*1e3, linetype=2) + geom_vline(xintercept=16550*1e3, linetype=2) + 
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 20924180)) + scale_y_continuous(labels=function(x)x*1e2, limits=c(0,1.05)) +
  facet_grid(~chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

X <- syn_windows %>% filter(chr == "X") %>% 
  ggplot(aes(x=mid, y=prop)) + 
  geom_point() + 
  geom_smooth(color="black", fill="light grey", span=0.75) +
  geom_vline(xintercept=6137*1e3, linetype=2) + geom_vline(xintercept=12480*1e3, linetype=2) +
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 17718942)) + scale_y_continuous(labels=function(x)x*1e2, limits=c(0,1.05)) +
  facet_grid(~chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

PlotList1 <- list(I, II, III, IV, V, X)

# CELEG and CBRIG gene positions 
synteny_table <- read.table("CELEG_CBRIG.txt", sep="\t", col.names=c("sp1", "sp1seq", "sp1chr", "sp1start", "sp1end", "sp1order", "sp2", "sp2seq", "sp2chr", "sp2start", "sp2end", "sp2order"))

I <- synteny_table %>% filter(sp1chr == sp2chr) %>% filter(sp1chr == "I") %>% 
  mutate(sp1pos = ((sp1start+sp1end)/2)) %>% 
  mutate(sp2pos = ((sp2start+sp2end)/2)) %>%
  ggplot(., aes(x=sp1pos, y=sp2pos)) + 
  geom_point(alpha=0.4, size=1, shape=16) + 
  geom_vline(xintercept=3858*1e3, linetype=2) + geom_vline(xintercept=11040*1e3, linetype=2) +
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 15072434)) + scale_y_continuous(labels=function(x)x/1e6) +  
  facet_wrap(~sp1chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

II <- synteny_table %>% filter(sp1chr == sp2chr) %>% filter(sp1chr == "II") %>% 
  mutate(sp1pos = ((sp1start+sp1end)/2)) %>% 
  mutate(sp2pos = ((sp2start+sp2end)/2)) %>%
  ggplot(., aes(x=sp1pos, y=sp2pos)) + 
  geom_point(alpha=0.4, size=1, shape=16) + 
  geom_vline(xintercept=4879*1e3, linetype=2) + geom_vline(xintercept=12020*1e3, linetype=2) +
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 15279421)) + scale_y_continuous(labels=function(x)x/1e6) +  
  facet_wrap(~sp1chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

III <- synteny_table %>% filter(sp1chr == sp2chr) %>% filter(sp1chr == "III") %>% 
  mutate(sp1pos = ((sp1start+sp1end)/2)) %>% 
  mutate(sp2pos = ((sp2start+sp2end)/2)) %>%
  ggplot(., aes(x=sp1pos, y=sp2pos)) + 
  geom_point(alpha=0.4, size=1, shape=16) + 
  geom_vline(xintercept=3722*1e3, linetype=2) + geom_vline(xintercept=10340*1e3, linetype=2) +
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 13783801)) + scale_y_continuous(labels=function(x)x/1e6) +  
  facet_wrap(~sp1chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

IV <- synteny_table %>% filter(sp1chr == sp2chr) %>% filter(sp1chr == "IV") %>% 
  mutate(sp1pos = ((sp1start+sp1end)/2)) %>% 
  mutate(sp2pos = ((sp2start+sp2end)/2)) %>%
  ggplot(., aes(x=sp1pos, y=sp2pos)) + 
  geom_point(alpha=0.4, size=1, shape=16) + 
  geom_vline(xintercept=3896*1e3, linetype=2) + geom_vline(xintercept=12970*1e3, linetype=2) +
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 17493829)) + scale_y_continuous(labels=function(x)x/1e6) +  
  facet_wrap(~sp1chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

V <- synteny_table %>% filter(sp1chr == sp2chr) %>% filter(sp1chr == "V") %>% 
  mutate(sp1pos = ((sp1start+sp1end)/2)) %>% 
  mutate(sp2pos = ((sp2start+sp2end)/2)) %>%
  ggplot(., aes(x=sp1pos, y=sp2pos)) + 
  geom_point(alpha=0.4, size=1, shape=16) + 
  geom_vline(xintercept=5897*1e3, linetype=2) + geom_vline(xintercept=16550*1e3, linetype=2) +
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 20924180)) + scale_y_continuous(labels=function(x)x/1e6) +  
  facet_wrap(~sp1chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

X <- synteny_table %>% filter(sp1chr == sp2chr) %>% filter(sp1chr == "X") %>% 
  mutate(sp1pos = ((sp1start+sp1end)/2)) %>% 
  mutate(sp2pos = ((sp2start+sp2end)/2)) %>%
  ggplot(., aes(x=sp1pos, y=sp2pos)) + 
  geom_point(alpha=0.4, size=1, shape=16) + 
  geom_vline(xintercept=6137*1e3, linetype=2) + geom_vline(xintercept=12480*1e3, linetype=2) +
  scale_x_continuous(labels=function(x)x/1e6, limits=c(0, 17718942)) + scale_y_continuous(labels=function(x)x/1e6) +  
  facet_wrap(~sp1chr) + theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0)

PlotList2 <- list(I, II, III, IV, V, X)

p1 <- PlotList1 %>% append(PlotList2) %>%
  ArrangeCowplot(.) + 
  plot_layout(nrow = 2)

## CINOP/CNIGO
## CINOP and CNIGO synteny windows
syn_windows <- read.table("CINOP_CNIGO.syn_windows", sep="\t", col.names=c("chr", "start", "end", "syn", "break", "prop"))
syn_windows <- syn_windows %>% mutate(mid = (start + end)/2)

PlotList1 <- 
  syn_windows$chr %>% 
  unique %>% 
  sort %>% 
  map(~syn_windows %>%
        filter(chr == .x) %>% 
        ggplot(aes(x=mid, y=prop)) + 
        geom_point() + 
        geom_smooth(color="black", fill="light grey", span=0.75) +
        scale_x_continuous(labels=function(x)x/1e6) +
        scale_y_continuous(labels=function(x)x*1e2, limits=c(0,1.05)) +
        facet_grid(~chr) +
        theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0))


## CINOP and CNIGO gene positions 
synteny_table <- read.table("CINOP_CNIGO.txt", sep="\t", col.names=c("sp1", "sp1seq", "sp1chr", "sp1start", "sp1end", "sp1order", "sp2", "sp2seq", "sp2chr", "sp2start", "sp2end", "sp2order"))

PlotList2 <- 
  synteny_table$sp1chr %>% 
  unique %>% 
  sort %>% 
  map(~synteny_table %>%
        filter(sp1chr == sp2chr) %>%
        filter(sp1chr == .x) %>% 
        mutate(sp1pos = ((sp1start+sp1end)/2)) %>% 
        mutate(sp2pos = ((sp2start+sp2end)/2)) %>%
        ggplot(., aes(x=sp1pos, y=sp2pos)) + 
        geom_point(alpha=0.4, size=1, shape=16) + 
        scale_x_continuous(labels=function(x)x/1e6) + 
        scale_y_continuous(labels=function(x)x/1e6) +  
        facet_wrap(~sp1chr) +         theme_bw() + theme(axis.title = element_blank()) + expand_limits(x=0))

p2 <- PlotList1 %>% append(PlotList2) %>%
  ArrangeCowplot(.) + 
  plot_layout(nrow = 2)


p <- p1/p2

ggsave("all_synteny-arms_centres.pdf", p, width=16, height=12)
ggsave("all_synteny-arms_centres.png", p, width=16, height=12)
