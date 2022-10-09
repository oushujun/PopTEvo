## PopTEvo R code for pan-genome TE analyses
## Shujun Ou (shujun.ou.1@gmail.com)
## 09/26/2022

#update.packages()
library(tidyr)
library(ggplot2)
library(wesanderson)
library(scales)
library(reshape2)
library(plyr)
library(grid)
library(gridExtra)
library(lattice)
#install.packages("devtools")
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(RColorBrewer)
library(vioplot)
library(ggridges)
library("stringr")
library(dplyr) # require dplyr 1.0.5 or lower by ggtree
library(ggsignif)
library(MKinfer)
#devtools::install_github("solatar/ggbrace")
library(ggbrace)
library(ggtext)
library(ggforce)
if (!require(coin)) {
  install.packages("coin")
}
if (!require(multcompView)) {
  install.packages("multcompView")
}
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
require(rcompanion)
library("ggVennDiagram")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(rvcheck)


datadir = "/Users/oushujun/My Drive/study/Maize research/NAM/TE/panArabidopsis/"
setwd(datadir)

########################
#### basic settings ####
########################
TE_colors = c('gray', rev(brewer.pal(11,'RdYlBu')[1:5]), 
              brewer.pal(11,'RdYlBu')[6], 
              rev(brewer.pal(11,'RdBu')[7:10]))

show_col(TE_colors)

TE_colors2 = c("LTR/Copia"="#D1E5F0", "LTR/Ty3"="#92C5DE", "LTR/unknown"="#4393C3", 
               "LINE/L1"="#2166AC",  "LINE/RTE"="#2166AC", "LINE/unknown"="#2166AC",
               "DNA/Helitron"="#FFFFBF", "TIR/Tc1_Mariner"="#A50026", "TIR/Mutator"="#D73027", 
               "TIR/PIF_Harbinger"="#F46D43", "TIR/CACTA"="#FDAE61", "TIR/hAT"="#FEE090", 
               "centromeric_repeat"="gray", "rDNA_spacer"="gray", "subtelomere"="gray", 
               "knob"="gray", "low_complexity"="gray")
show_col(TE_colors2)

TE_colors3 = c("LTR/Copia"="#92C5DE", "LTR/Ty3"="#92C5DE", "LTR/unknown"="#92C5DE", 
               "LINE/L1"="#2166AC",  "LINE/RTE"="#2166AC", "LINE/unknown"="#2166AC", "LINE"="#2166AC", 
               "DNA/Helitron"="#D73027", "TIR/Tc1_Mariner"="#FDAE61", "TIR/Mutator"="#FDAE61", 
               "TIR/PIF_Harbinger"="#FDAE61", "TIR/CACTA"="#FDAE61", "TIR/hAT"="#FDAE61", 
               "centromeric_repeat"="gray", "rDNA_spacer"="gray", "subtelomere"="gray", 
               "knob"="gray", "low_complexity"="gray")
show_col(TE_colors3)

rename_supfam = list("Gypsy_LTR_retrotransposon"="LTR/Ty3", 
                     "LTR/Gypsy"="LTR/Ty3",
                     "LTR/CRM"="LTR/Ty3",
                     'Gypsy'="LTR/Ty3",
                     "Copia_LTR_retrotransposon"="LTR/Copia", 
                     'Copia'="LTR/Copia", 
                     "LTR_retrotransposon"="LTR/unknown", 
                     "LTR/Retrovirus"="LTR/unknown",
                     "LTR/mixture"="LTR/unknown",
                     "LTR/Bel-Pao"="LTR/unknown",
                     "LTR_unknown"="LTR/unknown",
                     "helitron"="DNA/Helitron", 
                     "Helitron/unknown"="DNA/Helitron",
                     "PIF_Harbinger_TIR_transposon"="TIR/PIF_Harbinger",
                     "MITE/DTH"="TIR/PIF_Harbinger",
                     "DNA/DTH"="TIR/PIF_Harbinger",
                     "PIF_Harbinger"="TIR/PIF_Harbinger",
                     "hAT_TIR_transposon"="TIR/hAT", 
                     "DNA/DTA"="TIR/hAT", 
                     "MITE/DTA"="TIR/hAT",
                     'hAT'="TIR/hAT",
                     "Mutator_TIR_transposon"="TIR/Mutator", 
                     "TIR/MuDR_Mutator"="TIR/Mutator", 
                     "DNA/DTM"="TIR/Mutator",
                     "MITE/DTM"="TIR/Mutator",
                     "Mutator"="TIR/Mutator",
                     "CACTA_TIR_transposon"="TIR/CACTA", 
                     "TIR/EnSpm_CACTA"="TIR/CACTA",
                     "CACTA"="TIR/CACTA",
                     "DNA/DTC"="TIR/CACTA", 
                     "MITE/DTC"="TIR/CACTA",  
                     "Tc1_Mariner_TIR_transposon"="TIR/Tc1_Mariner", 
                     "MITE/DTT"="TIR/Tc1_Mariner", 
                     "DNA/DTT"="TIR/Tc1_Mariner", 
                     "Tc1_Mariner"="TIR/Tc1_Mariner", 
                     "TIR/Sola1"="TIR/unknown",
                     "TIR/Sola2"="TIR/unknown",
                     "TIR/Merlin"="TIR/unknown",
                     "TIR/Kolobok"="TIR/unknown",
                     "TIR/Novosib"="TIR/unknown",
                     "L1_LINE_retrotransposon"="LINE/L1", 
                     "RTE_LINE_retrotransposon"="LINE/RTE", 
                     "LINE_element"="LINE/unknown",
                     "Maverick/unknown"='unknown',
                     "Penelope/unknown"='unknown',
                     "mixture/mixture"='unknown',
                     "DIRS/unknown"='unknown',
                     "pararetrovirus/unknown"='unknown',
                     "CentC"="centromeric_repeat", 
                     "rDNA_intergenic_spacer_element"="rDNA/spacer", 
                     "subtelomere"="subtelomere", 
                     "knob"="knob", 
                     "low_complexity"="low_complexity",
                     'tRNA_SINE_retrotransposon'='SINE/tRNA')

####################################
#### Genome-level TE size stats ####
####################################
# import species
genome_info = read.csv('ncbi_dataset.50.csv')
genome_info$Assembly.Accession = sub('\\.[0-9]+', '', genome_info$Assembly.Accession)
LAI = read.table('Arabidopsis.panEDTA.LAI', header = F)
names(LAI) = c('genome', 'start', 'size', 'intact', 'total', 'rawLAI', 'LAI')

## percent TE in assembly
repeats = read.table('Arabidopsis.panEDTA.TE.anno.pcnt.txt', header=T, fill = T)
repeats$Ty3 = repeats$Gypsy #convert supfam name

# remove % and convert to numeric; remove _pcnt string in genome names
repeats[,-1] = apply(repeats[,-1], 2, function(x){as.numeric(sub("%", "", x, fixed=TRUE))})
repeats$Genome = sub("_pcnt", "", repeats$Genome, fixed=TRUE)
str(repeats)

# add more categories
repeats$DNA = repeats$helitron + repeats$CACTA + repeats$Mutator + repeats$PIF_Harbinger +
  repeats$Tc1_Mariner + repeats$hAT + repeats$DNA_transposon
repeats$nonLTR = repeats$L1_LINE + repeats$LINE_element + repeats$SINE_element + repeats$tRNA_SINE
repeats$LTR = repeats$Copia + repeats$Ty3 + repeats$LTR_unknown
repeats$nonTE = repeats$satellite_DNA + repeats$low_complexity + repeats$repeat_region

# gather plotting variables and change the variable plotting order
repeats_final = subset(repeats, !is.na(Total)) %>% 
  gather(variable, value, hAT, CACTA, PIF_Harbinger, Mutator, Tc1_Mariner,  
         helitron, nonLTR, Copia, Ty3, LTR_unknown, nonTE)
repeats_final$variable = as.factor(repeats_final$variable) #convert to character to factor
repeats_final$variable = factor(repeats_final$variable, 
                                levels = c("nonTE", "hAT", "CACTA", "PIF_Harbinger", "Mutator", 
                                           "Tc1_Mariner", "helitron", "nonLTR", 
                                           "LTR_unknown", "Ty3", "Copia"))
repeats_final$species = as.factor(genome_info$Organism.Name[match(repeats_final$Genome, genome_info$Assembly.Accession)])

# reorder plotting order
repeats_final$Genome = factor(repeats_final$Genome, levels = genome_info$Assembly.Accession)
repeats_final$LAI = LAI$LAI[match(repeats_final$Genome, LAI$genome)]
repeats_pcnt_plot = repeats_final %>%
  ggplot() + 
  geom_bar(aes(x = Genome, y = value, fill = variable), stat = "identity",colour="black") +
  geom_point(aes(x = Genome, y = LAI)) +
  scale_fill_manual(values=TE_colors) +
  scale_y_continuous(breaks = seq(0, 35, by = 5), 
                     sec.axis = sec_axis(~ ., name="LTR Assembly Index", breaks = seq(0, 35, by = 5))) +
  geom_point(aes(y = -1, x = Genome, shape = species), size = 2) +
  scale_shape_manual(values=0:nlevels(repeats_final$Genome)) +
  labs(x =" ", y = "Repeat content (%)",size=12, face="bold") +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95),
        legend.text = element_text(size=12),
        legend.title=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top", legend.box="vertical", 
        legend.margin=margin(), legend.text = element_text(size=10),
        plot.margin = margin(10, 10, 10, 25))
repeats_pcnt_plot

## total TE bp in assembly
repeats_bp = read.table('Arabidopsis.panEDTA.TE.anno.bp.txt', header=T, fill = T)
repeats_bp$Genome = sub("_bp", "", repeats_bp$Genome, fixed=TRUE)
repeats_bp = repeats_bp[!names(repeats_bp) %in% c('centromeric_repeat', 'knob', 'low_complexity', 'rDNA_spacer', 
                                                  'satellite_DNA', 'repeat_region', 'subtelomere', 'Total')]
repeats_bp$TE = rowSums(repeats_bp[,-1])

## total intact TE bp in assembly
intact = read.table('Arabidopsis.panEDTA.TE.intact.sum.txt', header=T)
intact$intact = intact$LTR + intact$TIR + intact$Helitron
intact$fragmented = repeats_bp$TE - intact$intact
intact$total = intact$intact + intact$fragmented
str(intact)
intact$Genome
rm(repeats_bp)

#summarize intact and fragmented percentage
intact_pcntTE = data.frame(intact[,-1]*100/(intact$total), Genome=intact$Genome) #intact/TE, frag/TE
summary(intact_pcntTE$intact)
summary(intact_pcntTE$fragmented)
std <- function(x) sd(x)/sqrt(length(x)) #standard error
sd(intact_pcntTE$intact)
sd(intact_pcntTE$fragmented)
std(intact_pcntTE$intact)
std(intact_pcntTE$fragmented)

# plot the %intact and %homo in %TE
intact_pcntTE_final = intact_pcntTE %>% gather(variable, value, intact, fragmented)
intact_pcntTE_final$variable = as.factor(intact_pcntTE_final$variable) #convert to character to factor
intact_pcntTE_final$variable = factor(intact_pcntTE_final$variable, 
                                      levels = c("intact", "fragmented"))

# plot the %intact and %frag to the total TE size
intact_pcntTE_plot = intact_pcntTE_final %>%
  ggplot(aes(x = Genome, y = value, fill = variable)) + 
  geom_bar(stat = "identity", color="black") +
  scale_fill_manual(values=c("#92C5DE", "gray")) +
  labs(x =" ", y = "TE length (%)",size=12, face="bold") +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95),
        legend.text = element_text(size=12),
        legend.title=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="top")
intact_pcntTE_plot



## Plot Figure to file
pdf("Suppl Figure x.pdf", width=11,height=8,pointsize=12, paper='special')
repeats_pcnt_plot
dev.off()
