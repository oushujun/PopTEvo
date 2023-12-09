#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)

# READ INPUT
stat <- read.csv("../data/nlr_information.tsv",header = T, sep = "\t")
te <- read.csv("../data/ltr_neighborhood_age.tsv",header = T, sep = "\t")
temp <- read.csv("../data/ltr_neighborhood_temperate_tropical-amplifying.tsv",header = T, sep = "\t") 
trop <- read.csv("../data/ltr_neighborhood_tropical_tropical-amplifying.tsv",header = T, sep = "\t") 
dist <- read.csv("../data/distance_ltr_neighborhood.tsv",header = T, sep = "\t") 

# PROCESS INPUT DATA FOR PLOTTING

# Summary stats
stat$size <- stat$end - stat$start
allchr <- stat
allchr$chr <- "all"
statall <- rbind(stat,allchr)
statall $chr <- factor(statall $chr, levels=c("all", "chr1", "chr2","chr3", "chr4", "chr5", "chr6","chr7", "chr8","chr9","chr10"))
statnonnlr <- statall[statall$gene_group!="nlr",]
statall <- statall[statall$gene_group=="nlr",]

chrstats <- statall %>% 
  group_by(chr,nam_line,cluster) %>% 
  summarise(sum_genecount = sum(genecount))

cluststats <- statall %>% 
  group_by(chr,nam_line,cluster) %>% 
  summarise(clustercount = n())

# Sum up all of the lines
tesum <- te %>% group_by(overlap,test,k) %>%
  summarize(Sum_bgd_young = sum(bgd_young),
            Sum_bgd_old = sum(bgd_old),
            Sum_bgd_notyoung = sum(bgd_notyoung),
            Sum_fgd_young = sum(fgd_young),
            Sum_fgd_old = sum(fgd_old),
            Sum_fgd_notyoung = sum(fgd_notyoung),
            young_enrichment =(sum(fgd_young)/sum(fgd_old))/ (sum(bgd_young)/sum(bgd_old)),
            sd_enrichment = sd((fgd_young/fgd_old)/(bgd_young/bgd_old)),
            young_enrichment_max = max((fgd_young/fgd_old)/(bgd_young/bgd_old)),
            young_enrichment_min = min((fgd_young/fgd_old)/(bgd_young/bgd_old))
  )

tesum_main <- tesum[tesum$test=="tNLR100kb_tnonNLR100kb" | tesum$test=="tNLR100kb_snonNLR"| tesum$test=="tnonNLR100kb_snonNLR" ,]

# Analysing tropical amplification

ampnam <- merge(temp,trop,by = c("overlap","test","k"))

ampsum <- ampnam %>% group_by(overlap,test,k) %>%
  summarize(Sum_temp_bgd_notamp = sum(temp_bgd_notamp),
            Sum_temp_bgd_amp = sum(temp_bgd_amp),
            Sum_temp_fgd_notamp = sum(temp_fgd_notamp),
            Sum_temp_fgd_amp = sum(temp_fgd_amp),
            Sum_trop_bgd_notamp = sum(trop_bgd_notamp),
            Sum_trop_bgd_amp = sum(trop_bgd_amp),
            Sum_trop_fgd_notamp = sum(trop_fgd_notamp),
            Sum_trop_fgd_amp = sum(trop_fgd_amp),
            sd_enrichment = sd((trop_fgd_amp/trop_bgd_amp)/(temp_fgd_amp/temp_bgd_amp))
            )

ampsumsub <- ampsum[, c("overlap", "test","k","sd_enrichment")]

tempsum <- temp %>% group_by(overlap,test,k) %>%
  summarize(Sum_temp_bgd_notamp = sum(temp_bgd_notamp),
            Sum_temp_bgd_amp = sum(temp_bgd_amp),
            Sum_temp_fgd_notamp = sum(temp_fgd_notamp),
            Sum_temp_fgd_amp = sum(temp_fgd_amp)
  )
tropsum <- trop %>% group_by(overlap,test,k) %>%
  summarize(Sum_trop_bgd_notamp = sum(trop_bgd_notamp),
            Sum_trop_bgd_amp = sum(trop_bgd_amp),
            Sum_trop_fgd_notamp = sum(trop_fgd_notamp),
            Sum_trop_fgd_amp = sum(trop_fgd_amp)
  )

amp <- merge(tempsum,tropsum,by = c("overlap","test","k"))

ampmut <- amp %>% group_by(overlap,test,k) %>%
  mutate(amp_enrichment =(Sum_trop_fgd_amp/Sum_trop_bgd_amp)/(Sum_temp_fgd_amp/Sum_temp_bgd_amp)
  )

ampmutsd <- merge(ampmut,ampsumsub,by = c("overlap","test","k"))
ampmutsd_main <- ampmutsd[ampmutsd$test=="tNLR100kb_tnonNLR100kb" | ampmutsd$test=="tNLR100kb_snonNLR"| ampmutsd$test=="tnonNLR100kb_snonNLR" ,]

# Add distance
dist <- dist[dist$clust_group!="50kb" & dist$info=="OVERLAPDIST",]
# Sum up all of the lines
distsum <- dist %>% group_by(k) %>%
  summarize(mean_dist = mean(dist_bp))

# PLOT FIGURES

# Figure 1
# Summary statistics for the NLR annotation and clustering in NAM lines.

p1 <- ggplot(statall[statall$cluster=="tandem_100kb",],aes(chr,size/1000,fill="#7570b3"))+
  geom_boxplot() + theme_bw() + theme(text = element_text(size=20),legend.position="none")+
  ylab("NLR cluster size (kb)")+
  xlab("")+
  scale_fill_manual(values = c("#7570b3"))

p2 <- ggplot(statall[statall$cluster=="tandem_100kb",],aes(chr,genecount,fill="#7570b3"))+
  geom_boxplot() + theme_bw() + theme(text = element_text(size=20),legend.position="none")+
  ylab("NLR cluster gene count")+
  xlab("")+
  scale_fill_manual(values = c("#7570b3"))

p3 <- ggplot(chrstats[chrstats$cluster!="tandem_50kb",],aes(chr,sum_genecount,fill=cluster))+
  geom_boxplot() + theme_bw()+ theme(text = element_text(size=20),legend.position="none")+
  ylab("NLR gene count per NAM line")+
  xlab("")+
  scale_fill_manual(values=c("#99d8c9","#2ca25f"))

p4 <- ggplot(cluststats[cluststats$cluster=="tandem_100kb",],aes(chr,clustercount,fill="#7570b3"))+
  geom_boxplot() + theme_bw() + theme(text = element_text(size=20),legend.position="none")+
  ylab("NLR clusters per NAM line")+
  xlab("")+
  scale_fill_manual(values = c("#7570b3"))

legend_b <- get_legend(
  p3 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(p1,p2,p3,p4)
ggsave("Figure1.pdf", 
       width = 36,
       height = 18,
       units = c("cm"),
       dpi = 300)

# Figure 2
# Fold enrichment of Young LTR families compared to Old LTR families is increased in the neighborhood of NLR clusters compared to non-NLR gene clusters.
ggplot(tesum_main[tesum_main$overlap=="OVERLAP",],aes(x=k,y=young_enrichment,color=test,fill=test)) + geom_line()+
  geom_ribbon(aes(ymin=young_enrichment-sd_enrichment,ymax=young_enrichment+sd_enrichment),alpha=0.2)+
  theme_bw()+
  xlim(5,50) +
  scale_x_continuous("TEs included in neighborhood (k)", 
                     sec.axis = sec_axis(~ . * 22226.43 / 1000,name = "Maximum TE distance (kbp)"),
                     limits = c(5,50),
                     breaks = seq(5, 50, 5),
                     expand = c(0.01, 0.01)
  ) +
  ylab("Fold enrichment of young TEs")+
  xlab("TEs included in neighborhood (k)")+
  theme(text = element_text(size= 25),legend.title= element_blank())+
  scale_fill_manual(values = c("#2ca25f","#99d8c9","darkgrey")) +
  scale_color_manual(values = c("#2ca25f","#99d8c9","darkgrey"))

ggsave("Figure2.pdf", 
       width = 36,
       height = 18,
       units = c("cm"),
       dpi = 300)


# Figure 3
# Young LTR families are not enriched in the neighborhood of NLR singleton genes compared to non-NLR gene clusters.

tesum_supp1 <- tesum[tesum$test=="sNLR_snonNLR" | tesum$test=="tNLR100kb_sNLR" ,]
ggplot(tesum_supp1[tesum_supp1$overlap=="OVERLAP",],aes(x=k,y=young_enrichment,color=test,fill=test)) + geom_line()+
  geom_ribbon(aes(ymin=young_enrichment-sd_enrichment,ymax=young_enrichment+sd_enrichment),alpha=0.2)+
  theme_bw()+
  xlim(5,50) +
  scale_x_continuous("TEs included in neighborhood (k)", 
                     sec.axis = sec_axis(~ . * 22226.43 / 1000,name = "Maximum TE distance (kbp)"),
                     limits = c(5,50),
                     breaks = seq(5, 50, 5),
                     expand = c(0.01, 0.01)
  ) +
  ylab("Fold enrichment of young TEs")+
  xlab("TEs included in neighborhood (k)")+
  theme(text = element_text(size=25),legend.title= element_blank())+
  scale_fill_manual(values = c("darkgrey","#2ca25f")) +
  scale_color_manual(values = c("darkgrey","#2ca25f"))

ggsave("Figure3.pdf", 
       width = 28,
       height = 18,
       units = c("cm"),
       dpi = 300)

# Figure 4
# Young LTR families are enriched in NLR clusters at a clustering distance threshold of 50 kb

tesum_supp2 <- tesum[tesum$test=="tnonNLR50kb_snonNLR" | tesum$test=="tNLR50kb_sNLR" | tesum$test=="tNLR50kb_snonNLR"| tesum$test=="tNLR50kb_tnonNLR50kb" ,]
ggplot(tesum_supp2[tesum_supp2$overlap=="OVERLAP",],aes(x=k,y=young_enrichment,color=test,fill=test)) + geom_line()+
  geom_ribbon(aes(ymin=young_enrichment-sd_enrichment,ymax=young_enrichment+sd_enrichment),alpha=0.2)+
  theme_bw()+
  xlim(5,50) +
  scale_x_continuous("TEs included in neighborhood (k)", 
                     sec.axis = sec_axis(~ . * 22226.43 / 1000,name = "Maximum TE distance (kbp)"),
                     limits = c(5,50),
                     breaks = seq(5, 50, 5),
                     expand = c(0.01, 0.01)
  ) +
  ylab("Fold enrichment of young TEs")+
  xlab("TEs included in neighborhood (k)")+
  theme(text = element_text(size= 25),legend.title= element_blank())+
  scale_fill_manual(values = c("#7570b3","#2ca25f","#99d8c9","darkgrey")) +
  scale_color_manual(values = c("#7570b3","#2ca25f","#99d8c9","darkgrey"))

ggsave("Figure4.pdf", 
       width = 28,
       height = 18,
       units = c("cm"),
       dpi = 300)

# Figure 5
# Tropical amplifying LTR families in the neighborhood of NLR clusters are not significantly enriched in tropical lines compared to temperate lines.
ggplot(ampmutsd_main[ampmutsd_main$overlap=="OVERLAP",],aes(x=k,y=amp_enrichment,color=test,fill=test)) + geom_line()+
  geom_ribbon(aes(ymin=amp_enrichment-sd_enrichment,ymax=amp_enrichment+sd_enrichment),alpha=0.2)+
  theme_bw()+
  xlim(5,50) +
  scale_x_continuous("TEs included in neighborhood (k)", 
                     sec.axis = sec_axis(~ . * 22226.43 / 1000,name = "Maximum TE distance (kbp)"),
                     limits = c(5,50),
                     breaks = seq(5, 50, 5),
                     expand = c(0.01, 0.01)
  ) +
  ylab("Fold enrichment of tropical-amplified TEs in Tropicals")+
  xlab("TEs included in neighborhood (k)")+
  theme(text = element_text(size= 25),legend.title= element_blank())+
  scale_fill_manual(values = c("#2ca25f","#99d8c9","darkgrey")) +
  scale_color_manual(values = c("#2ca25f","#99d8c9","darkgrey"))

ggsave("Figure5.pdf", 
       width = 28,
       height = 24,
       units = c("cm"),
       dpi = 300)

