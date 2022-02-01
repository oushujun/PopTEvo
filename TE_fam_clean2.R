## PopTEvo R code for pan-genome TE analyses
## Shujun Ou (shujun.ou.1@gmail.com)
## 02/01/2022

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
library(dplyr)
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
#if(!require(ggtree)){BiocManager::install("ggtree", force = TRUE)}

datadir = "/Users/oushujun/Google\ Drive/My\ Drive/study/Maize\ research/NAM/TE"
setwd(datadir)
load("./TE_fam_PCA2.RData")

########################
#### basic settings ####
########################

NAM_ID = c("B73", "B97", "Ky21", "M162W", "Ms71", "Oh43", "Oh7B", "M37W", "Mo18W", "Tx303",
           "HP301", "P39", "Il14H", "CML52", "CML69", "CML103", "CML228", "CML247", "CML277", 
           "CML322", "CML333", "Ki3", "Ki11", "NC350", "NC358", "Tzi8")

NAM_temp = c("B73", "B97", "Ky21", "M162W", "Ms71", "Oh43", "Oh7B")
NAM_flnt = c("HP301", "P39", "Il14H")
NAM_adx = c("M37W", "Mo18W", "Tx303")
NAM_nontrop = c(NAM_temp, NAM_flnt)
NAM_trop = c("CML52", "CML69", "CML103", "CML228", "CML247", "CML277", "CML322", "CML333", 
             "Ki3", "Ki11", "NC350", "NC358", "Tzi8")
NAM = data.frame(genome = NAM_ID, group = c(rep('temp', 7), rep('admx', 3), rep('temp', 3), rep('trop', 13)))

NAM_colors = c("B73"="goldenrod1", "B97"="royalblue", "CML103"="limegreen", 
               "CML228"="limegreen", "CML247"="limegreen", "CML277"="limegreen", 
               "CML322"="limegreen", "CML333"="limegreen", "CML52"="limegreen", 
               "CML69"="limegreen", "HP301"="orchid", "Il14H"="orangered", 
               "Ki11"="limegreen", "Ki3"="limegreen", "Ky21"="royalblue", 
               "M162W"="royalblue", "M37W"="gray47", "Ms71"="royalblue", 
               "Mo18W"="gray47", "NC350"="limegreen", "NC358"="limegreen", 
               "Oh43"="royalblue", "Oh7B"="royalblue", "P39"="orangered", 
               "Tx303"="gray47", "Tzi8"="limegreen")
#reorder NAM_colors based on NAM_ID
NAM_colors = NAM_colors[order(match(names(NAM_colors), NAM_ID))]
TE_colors = c('gray', rev(brewer.pal(11,'RdYlBu')[1:5]), 
              brewer.pal(11,'RdYlBu')[6], 
              rev(brewer.pal(11,'RdBu')[7:10]))

show_col(TE_colors)

TE_colors2 = c("LTR/Copia"="#D1E5F0", "LTR/Gypsy"="#92C5DE", "LTR/unknown"="#4393C3", 
               "LINE/L1"="#2166AC",  "LINE/RTE"="#2166AC", "LINE/unknown"="#2166AC",
               "DNA/Helitron"="#FFFFBF", "TIR/Tc1_Mariner"="#A50026", "TIR/Mutator"="#D73027", 
               "TIR/PIF_Harbinger"="#F46D43", "TIR/CACTA"="#FDAE61", "TIR/hAT"="#FEE090", 
               "centromeric_repeat"="gray", "rDNA_spacer"="gray", "subtelomere"="gray", 
               "knob"="gray", "low_complexity"="gray")
show_col(TE_colors2)

TE_colors3 = c("LTR/Copia"="#92C5DE", "LTR/Gypsy"="#92C5DE", "LTR/unknown"="#92C5DE", 
               "LINE/L1"="#2166AC",  "LINE/RTE"="#2166AC", "LINE/unknown"="#2166AC", "LINE"="#2166AC", 
               "DNA/Helitron"="#D73027", "TIR/Tc1_Mariner"="#FDAE61", "TIR/Mutator"="#FDAE61", 
               "TIR/PIF_Harbinger"="#FDAE61", "TIR/CACTA"="#FDAE61", "TIR/hAT"="#FDAE61", 
               "centromeric_repeat"="gray", "rDNA_spacer"="gray", "subtelomere"="gray", 
               "knob"="gray", "low_complexity"="gray")
show_col(TE_colors3)

rename_supfam = list("Gypsy_LTR_retrotransposon"="LTR/Gypsy", 
                     "LTR/CRM"="LTR/Gypsy",
                     "Copia_LTR_retrotransposon"="LTR/Copia", 
                     "LTR_retrotransposon"="LTR/unknown", 
                     "LTR/Retrovirus"="LTR/unknown",
                     "LTR/mixture"="LTR/unknown",
                     "LTR/Bel-Pao"="LTR/unknown",
                     "helitron"="DNA/Helitron", 
                     "Helitron/unknown"="DNA/Helitron",
                     "PIF_Harbinger_TIR_transposon"="TIR/PIF_Harbinger",
                     "MITE/DTH"="TIR/PIF_Harbinger",
                     "DNA/DTH"="TIR/PIF_Harbinger",
                     "hAT_TIR_transposon"="TIR/hAT", 
                     "DNA/DTA"="TIR/hAT", 
                     "MITE/DTA"="TIR/hAT",
                     "Mutator_TIR_transposon"="TIR/Mutator", 
                     "TIR/MuDR_Mutator"="TIR/Mutator", 
                     "DNA/DTM"="TIR/Mutator",
                     "MITE/DTM"="TIR/Mutator",
                     "CACTA_TIR_transposon"="TIR/CACTA", 
                     "TIR/EnSpm_CACTA"="TIR/CACTA",
                     "DNA/DTC"="TIR/CACTA", 
                     "MITE/DTC"="TIR/CACTA",  
                     "Tc1_Mariner_TIR_transposon"="TIR/Tc1_Mariner", 
                     "MITE/DTT"="TIR/Tc1_Mariner", 
                     "DNA/DTT"="TIR/Tc1_Mariner", 
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
                     "low_complexity"="low_complexity")




######################################
#### Read and process family data ####
######################################
fam = read.table('NAM.196/NAM.EDTA1.9.6.MTEC02052020.TE.v1.1.anno.sum.fam.bp', header=T)
fam_int = read.table('NAM.196/NAM.EDTA1.9.6.MTEC02052020.TE.v1.intact.sum.fam.bp', header=T)
fam_homo = read.table('NAM.196/NAM.EDTA1.9.6.MTEC02052020.TE.v1.homo.sum.fam.bp', header=T)
class = read.table('NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa.list', header=F)
colnames(class) = c('TE', "Supfam") #add column header

# rename classes
table(class$Supfam)
class$Supfam = dplyr::recode(class$Supfam, !!!rename_supfam)
table(class$Supfam)

# check if TEs are classified uniquely
length(class$TE)
length(unique(class$TE))

# fix genome names
fam$TE_fam = sub("Oh7b", "Oh7B", fam$TE_fam)
names(fam)[names(fam) == 'Oh7b'] = 'Oh7B'
fam_int$TE_fam = sub("Oh7b", "Oh7B", fam_int$TE_fam)
names(fam_int)[names(fam_int) == 'Oh7b'] = 'Oh7B'
fam_homo$TE_fam = sub("Oh7b", "Oh7B", fam_homo$TE_fam)
names(fam_homo)[names(fam_homo) == 'Oh7b'] = 'Oh7B'
str(fam)

#remove AB10
drops <- c("B73Ab10","B73_AB10")
fam = fam[ , !(names(fam) %in% drops)]
fam_int = fam_int[ , !(names(fam_int) %in% drops)]
fam_homo = fam_homo[ , !(names(fam_homo) %in% drops)]
colnames(fam) #check if Ab10 is removed
colnames(fam_int)
sum(str_detect(fam_int$TE_fam, ":", negate = T))

# transpose, and add back col and row names
library(data.table)
tr_fam <- transpose(fam[-1]) # remove the first row, $TE_fam 
rownames(tr_fam) = colnames(fam[-1])
colnames(tr_fam) = as.matrix(transpose(fam[1]))[1,]

tr_fam_int <- transpose(fam_int[-1]) # remove the first row, $TE_fam 
rownames(tr_fam_int) = colnames(fam_int[-1])
colnames(tr_fam_int) = as.matrix(transpose(fam_int[1]))[1,]
# check size
rowSums(tr_fam_int[,colnames(tr_fam_int) %in% Fam_count$TE])
colSums(fam_int[fam_int$TE_fam %in% Fam_count$TE, -1])

tr_fam_homo <- transpose(fam_homo[-1]) # remove the first row, $TE_fam 
rownames(tr_fam_homo) = colnames(fam_homo[-1])
colnames(tr_fam_homo) = as.matrix(transpose(fam_homo[1]))[1,]

#str(tr_fam)
fam = tr_fam
fam = fam[,apply(fam, 2, mean) != 0] #remove all 0 rows
rm(tr_fam)
fam_int = tr_fam_int
fam_int = fam_int[,apply(fam_int, 2, mean) != 0] #remove all 0 rows
rm(tr_fam_int)
fam_homo = tr_fam_homo
fam_homo = fam_homo[,apply(fam_homo, 2, mean) != 0] #remove all 0 rows
rm(tr_fam_homo)

# convert bp to mb
fam = fam/1000000
old_fam = fam
fam_int = fam_int/1000000
fam_homo = fam_homo/1000000
rowSums(fam_int[,colnames(fam_int) %in% Fam_count$TE])

# checked fam_homo, fam_int, no split families
length(str_subset(names(fam_homo), "TE_"))
length(unique(sub("_LTR|_INT", '', str_subset(names(fam_homo), "TE_"))))
names(fam_homo) = sub("_LTR|_INT", '', names(fam_homo))

length(str_subset(names(fam_int), "TE_"))
length(unique(sub("_LTR|_INT", '', str_subset(names(fam_int), "TE_"))))
names(fam_int) = sub("_LTR|_INT", '', names(fam_int))
rowSums(fam_int[,colnames(fam_int) %in% Fam_count$TE])

# combine LTR families split in two columns, this step take like 20 mins..
library(stringr)
length(str_subset(names(fam), "TE_"))
length(unique(sub("_LTR|_INT", '', str_subset(names(fam), "TE_"))))
fam_names = str_subset(names(fam), "TE_")
fam_all = sub("_LTR|_INT", '', str_subset(names(fam), "TE_"))
fam_dup = names(table(fam_all)[table(fam_all)>1])
length(fam_dup)
remove_list = c()
for (i in 1:length(fam_dup)){
  this_id = str_subset(fam_names, fam_dup[i])
  remove_list = c(remove_list, this_id)
  fam$this_fam = rowSums(fam[this_id]) #combine cols
  fam = fam[, !names(fam) %in% this_id] #remove cols
  names(fam)[names(fam) == 'this_fam'] = fam_dup[i] #rename new col
}
names(fam) = sub("_LTR|_INT", '', names(fam)) #convert the rest of names

# convert bp to TE space percentage
fam_pcnt = fam/rowSums(fam)


############################
#### Pan TE family stat ####
############################
## Remove entries that are not in the pan-TE library
fam_rare = colSums(fam[str_detect(colnames(fam), ":")])
length(fam_rare)
sum(fam_rare)/nrow(fam) #rare TE avg NAM, 39.43013
sum(colSums(fam))/length(rownames(fam)) #total TE avg NAM

fam = fam[-grep(":", colnames(fam),)] # remove rare TEs
sum(colSums(fam))/length(rownames(fam)) #common TE avg NAM #1795.437
sum(fam)/26 #common TE avg NAM #1795.437
dim(fam)

## count family occurrance
str(fam)
levels=c("0")
non_0_count <- 26 - sapply(levels, function(x)colSums(fam=="0")) #count occurrences of 0 in each column
non_0_count_int <- 26 - sapply(levels, function(x)colSums(fam_int=="0")) #count occurrences of 0 in each column
non_0_count_homo <- 26 - sapply(levels, function(x)colSums(fam_homo=="0")) #count occurrences of 0 in each column

# count and size stats
length(str_subset(colnames(fam), ":", negate = T))
length(str_subset(colnames(fam), ":", negate = F))
length(colnames(fam))
sum(colSums(fam)/non_0_count<0.1) #a list of TE fam <100kb on avg
sum(colSums(fam)/non_0_count>=0.1) #a list of TE fam >100kb on avg
sum(colSums(fam)/non_0_count>=0.1)/27228
sum(colSums(fam))/26 # total size of common fams on avg
total_TE = 1834.867 #rare + common avg NAM
sum(rowSums(fam[,colSums(fam)/non_0_count<0.1]))/26 #total size of fams < 0.1M, 104.0061Mb
sum(colSums(fam)[colSums(fam)/non_0_count>=0.1]/26) #NAM avg size of TE fam >100kb, 1691.431Mb
rowSums(fam[,colSums(fam)/non_0_count>=0.1])/total_TE
std = function(x) sd(x)/sqrt(length(x))
std(rowSums(fam[,colSums(fam)/non_0_count>=0.1])/total_TE)
std(rowSums(fam[,colSums(fam)/non_0_count<0.1])/total_TE)
sd(rowSums(fam[,colSums(fam)/non_0_count>=0.1])/total_TE)
sd(rowSums(fam[,colSums(fam)/non_0_count<0.1])/total_TE)
length(intersect(unique(sub('_LTR|_INT', '', class[str_detect(class$Supfam, 'LTR'),]$TE)), names(fam))) #total number of LTR fams
sum(colSums(fam[intersect(unique(sub('_LTR|_INT', '', class[str_detect(class$Supfam, 'LTR'),]$TE)), names(fam))]))/26 #total size of LTR fams
sum(colSums(fam[,LTR_fam_info$id])/26) #NAM avg size of classifiable fam (1058 fams)
sum(colnames(fam) %in% sub('_LTR|_INT', '', class[str_detect(class$Supfam, 'LTR'), ]$TE) & !colnames(fam) %in% LTR_fam_info$id) #total number of fragmented fams
sum(colSums(fam[colnames(fam) %in% sub('_LTR|_INT', '', class[str_detect(class$Supfam, 'LTR'), ]$TE) & !colnames(fam) %in% LTR_fam_info$id]))/26 #total size of fragmented fams
dim(subset(LTR, id %in% as.character(subset(Fam_count, count<25)$TE))) #446 intact LTRs present in less than 26 genomes
dim(subset(LTR, id %in% as.character(subset(Fam_count, count>25)$TE))) #1.3 million lTRs present in all 26 genomes

## add family and superfamily info and mean family size info
Fam_count = data.frame(count=non_0_count, TE=colnames(fam),
                       Supfam=class[match(colnames(fam), sub('_INT|_LTR', '', class$TE)),]$Supfam,
                       Size=colSums(fam)/non_0_count)
colnames(Fam_count) = c('count', 'TE', 'Supfam', 'Size')
dim(Fam_count)

Fam_count_int = data.frame(count=non_0_count_int, TE=colnames(fam_int),
                       Supfam=class[match(colnames(fam_int), sub('_INT|_LTR', '', class$TE)),]$Supfam,
                       Size=colSums(fam_int)/non_0_count_int)
colnames(Fam_count_int) = c('count', 'TE', 'Supfam', 'Size')
dim(Fam_count_int)

Fam_count_homo = data.frame(count=non_0_count_homo, TE=colnames(fam_homo),
                       Supfam=class[match(colnames(fam_homo), sub('_INT|_LTR', '', class$TE)),]$Supfam,
                       Size=colSums(fam_homo)/non_0_count_homo)
colnames(Fam_count_homo) = c('count', 'TE', 'Supfam', 'Size')
dim(Fam_count_homo)

## filter out nonTE categories
Fam_count = Fam_count[which(Fam_count$Supfam!="centromeric_repeat" & Fam_count$Supfam!="knob" &
                              Fam_count$Supfam!="low_complexity" & Fam_count$Supfam!="rDNA_intergenic_spacer_element" & 
                              Fam_count$Supfam!="subtelomere"),]
Fam_count$Supfam = factor(Fam_count$Supfam) #remove unused factor levels
levels(Fam_count$Supfam)
str(Fam_count)

## replace LINEs to nonLTR
nonLTR = c("LINE/L1", "LINE/RTE", "LINE/unknown")
Fam_count$Supfam <- as.character(Fam_count$Supfam)
Fam_count$Supfam <- factor(with(Fam_count, replace(Supfam, Supfam %in% nonLTR, "nonLTR")))
levels(Fam_count$Supfam)

# add family length info
fam_len = read.table('NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa.info', header = T, sep = '\t')
Fam_count$fam_len = fam_len$fam_len[match(Fam_count$TE, fam_len$id)]
rm(fam_len)

## Total number of TE families
length(Fam_count$count)

## Number of families that are present in all 26 NAM parents
Fam_count$TE = as.factor(Fam_count$TE)
sum = length(levels(Fam_count$TE[-grep(":", Fam_count$TE)])) #sum of TEs used in the NAM.EDTA.clean library
table(Fam_count$count)["26"]
table(Fam_count$count)["26"]/27228
table(Fam_count$count)["26"]/length(Fam_count$count)
table(Fam_count$count)["26"]/sum

## more than 20 times
sum(table(Fam_count$count)[20:26])
sum(table(Fam_count$count)[20:26])/27228
sum(table(Fam_count$count)[20:26])/length(Fam_count$count)
sum(table(Fam_count$count)[20:26])/sum

## Number of families that are only present in one NAM parent
table(Fam_count$count)["1"]
table(Fam_count$count)["1"]/27228
table(Fam_count$count)["1"]/length(Fam_count$count)
table(Fam_count$count)["1"]/sum
table(Fam_count[-grep(":", Fam_count$TE),]$count)["1"]/sum
table(Fam_count[-grep(":", Fam_count$TE),]$count)["1"]/27228

## other summaries
sum(table(Fam_count$count)[24:25])/27228
sum(table(Fam_count$count)[2:23])/27228
sum(table(Fam_count$count)[24:26])/27228

## histogram of freq
plotorder <- c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon", "LTR_retrotransposon", 
               "nonLTR", "helitron", "CACTA_TIR_transposon", "hAT_TIR_transposon", "Mutator_TIR_transposon", 
               "PIF_Harbinger_TIR_transposon", "Tc1_Mariner_TIR_transposon")
plotorder = dplyr::recode(plotorder, !!!rename_supfam)
Fam_count <- arrange(transform(Fam_count,
                               Supfam=factor(Supfam,levels=plotorder)),Supfam)

Fam_count$count = as.numeric(Fam_count$count)
Fam_count_freq_plot = ggplot(Fam_count, aes(count)) +
  geom_histogram(binwidth=1, center=0.5) +
  facet_wrap(~Supfam, scales = "free_y", ncol = 5) +
  labs(x ="Frequency (n = 26)", y = "Number of TE families") +
  theme(axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title=element_blank())
Fam_count_freq_plot

## readin bootstrap pan-TE data
pan_TE_bs = read.table('pan_TE_bootstrap1000.summary26.txt', header=F)
pan_TE_bs = data.frame(count=as.vector(t(pan_TE_bs)),Genome=rep(1:26, 1000))
pan_TE_bs$Genome = as.factor(pan_TE_bs$Genome)
length(pan_TE_bs)
mean(pan_TE_bs$count[pan_TE_bs$Genome==2])/21358
mean(pan_TE_bs$count[pan_TE_bs$Genome==5])/21358
mean(pan_TE_bs$count[pan_TE_bs$Genome==10])/21358
21358 - mean(pan_TE_bs$count[pan_TE_bs$Genome==10])

# add percentage to the second axis and change background colors
show_col(brewer.pal(8,'YlGn'))
color_levels = brewer.pal(8,'YlGn')
Total_flTE_fam = 21358
pan_TE_bs_plot_pcnt = ggplot(pan_TE_bs, aes(x=Genome, y=count)) + 
  annotation_raster(alpha(color_levels[2], 1), xmin = -Inf, xmax = Inf, 
                    ymin = -Inf, ymax = Total_flTE_fam*0.8) +
  annotation_raster(alpha(color_levels[3], 1), xmin = -Inf, xmax = Inf, 
                    ymin = Total_flTE_fam*0.8, ymax = Total_flTE_fam*0.9) +
  annotation_raster(alpha(color_levels[4], 1), xmin = -Inf, xmax = Inf, 
                    ymin = Total_flTE_fam*0.9, ymax = Total_flTE_fam*0.95) +
  annotation_raster(alpha(color_levels[5], 1), xmin = -Inf, xmax = Inf, 
                    ymin = Total_flTE_fam*0.95, ymax = Inf) +
  geom_violin() + #this layer at the top
  scale_x_discrete(breaks = c(1, 5, 10, 15, 20, 25)) +
  scale_y_continuous(sec.axis = sec_axis(~ . *100/Total_flTE_fam, breaks = c(75,80,85,90,95,100),name="% of TE families")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x ="Number of Genomes", y = "# of TE families") +
  theme(axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=12),
        legend.title=element_blank())
pan_TE_bs_plot_pcnt

## plot size of TE families
Fam_count$count = as.factor(Fam_count$count)
Fam_size_freq_int_plot_log = ggplot(Fam_count_int, aes(x=factor(count), y=Size)) + 
  geom_boxplot(outlier.size=0.8, width=0.7, fill="white") +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  labs(x ="Frequency (n = 26)", y = "TE family size (Mb)<br>Intact elements ")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.title.y = element_markdown(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title=element_blank())
Fam_size_freq_int_plot_log

Fam_size_freq_homo_plot_log = ggplot(Fam_count_homo, aes(x=factor(count), y=Size)) + 
  geom_boxplot(outlier.size=0.8, width=0.7, fill="white") +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  labs(x ="Frequency (n = 26)", y = "TE family size (Mb)<br>Fragmented elements")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.title.y = element_markdown(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title=element_blank())
Fam_size_freq_homo_plot_log

#number of families less than x Mb in total size
length(Fam_count[Fam_count$Size<0.01,]$Size)/27228 #10k
length(Fam_count[Fam_count$Size<0.1,]$Size)/27228 #100k
length(Fam_count[Fam_count$Size<1,]$Size)/27228 #1M

sum(Fam_count[Fam_count$Size<0.01,]$Size)/sum(Fam_count$Size)
sum(Fam_count[Fam_count$Size<0.1,]$Size)/sum(Fam_count$Size)
sum(Fam_count[Fam_count$Size<1,]$Size)/sum(Fam_count$Size)

#rare and large TE families
Fam_count$count = as.numeric(Fam_count$count)
Fam_count[Fam_count$Size>20 & Fam_count$count<=3,]
Fam_count[Fam_count$Size>0.02 & Fam_count$count<26,]
length(Fam_count[Fam_count$Size>0.1,]$count)

# plot rare LTR count data
rare_LTR2 = read.table('NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.LTR.rare.count2', header = F)
names(rare_LTR2) = c('genome', 'uniq', 'count', 'mean', 'sd')
rare_LTR2$genome = sub('\\..*', '', rare_LTR2$genome)
rare_LTR2 = rare_LTR2[rare_LTR2$genome != 'B73Ab10' & rare_LTR2$genome != 'B73_AB10', ] #rm Ab10
rare_LTR2$genome = sub("Oh7b", "Oh7B", rare_LTR2$genome)
rare_LTR2$genome = sub("IL14H", "Il14H", rare_LTR2$genome)
rare_LTR2$genome = sub("MS71", "Ms71", rare_LTR2$genome)
rare_LTR2$genome = factor(rare_LTR2$genome, NAM$genome)
unique(rare_LTR2$genome)
rare_LTR_plot2 = ggplot(rare_LTR2, aes(x=genome, y=mean)) + 
  geom_bar(stat = 'identity', alpha = 0.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position=position_dodge(.9)) +
  labs(x ="", y="Copy number of<br>unclassified LTR-RTs")+
  theme(axis.title.y = element_markdown(size=12, face="bold"),
        #  axis.title.y = element_text(),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10, angle=50,vjust=1,hjust=0.95)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
rare_LTR_plot2


############################
#### PCA ###################
############################

## PCA
names(fam)
rownames(fam)
fam$TE_00018010
dim(fam)
pca1 = prcomp(fam, scale. = F) #use unnormalized value to count in the effect of TE fam size
pca1$sdev       # sqrt of eigenvalues
pca1$rotation   # loadings
pca1$x          # PCs (aka scores)
scores = as.data.frame(pca1$x)
sum = summary(pca1)
str(sum)
PCvar_num = pca1$sdev^2/sum(pca1$sdev^2) #sum$importance, pcnt variations explained of PCs
PCvar = sprintf("%1.1f%%", 100*PCvar_num)
PCvar

# SNP PCA
NAM.SNP25k = as.data.frame(read.table('B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.homo.chr.gt.pres.25k.h', 
                        header = T) %>% select(all_of(NAM_ID)) %>% t(.))
snp_pca1 = prcomp(NAM.SNP25k, scale. = F) 
screeplot(snp_pca1)
snp_scores = as.data.frame(snp_pca1$x)
snp_sum = summary(snp_pca1)
snp_PCvar_num = snp_pca1$sdev^2/sum(snp_pca1$sdev^2) #sum$importance, pcnt variations explained of PCs
barplot(snp_PCvar_num, ylim=c(0,0.1), xlab='Principal components', ylab='Variations explained')
snp_PCvar = sprintf("%1.1f%%", 100*snp_PCvar_num)
snp_PCvar

#reorder columns based on NAM_colors
scores = scores[names(NAM_colors),]
rownames(scores)
snp_scores = snp_scores[names(NAM_colors),]
rownames(snp_scores)

# set ellipse group and plot
scores$group = rownames(scores)
scores[scores$group %in% tropical,]$group = 'trop'
scores[scores$group %in% non_tropical,]$group = 'temp'
scores[scores$group %in% NAM_adx,]$group = NA
PCA_biplot_group = ggplot(data=scores, aes(x=PC1, y=PC2)) +
  geom_hline(yintercept=0, colour="gray65") +
  geom_vline(xintercept=0, colour="gray65") +
  geom_mark_ellipse(data=subset(scores, !is.na(group)), expand = 0, aes(color=group), alpha = 0, 
                    linetype = 'dashed', size = 0.4, show.legend = F) +
  geom_text(aes(label=rownames(scores)),
            colour=NAM_colors, alpha=0.8, size=3) + 
  scale_color_manual(values=c('#53CACE', '#F5958F')) +
  xlab(paste("PC1 (", PCvar[1], ")", sep="")) + xlim(-9, 8) +
  ylab(paste("PC2 (", PCvar[2], ")", sep="")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
PCA_biplot_group

snp_scores$group = rownames(snp_scores)
snp_scores[snp_scores$group %in% tropical,]$group = 'trop'
snp_scores[snp_scores$group %in% non_tropical,]$group = 'temp'
snp_scores[snp_scores$group %in% NAM_adx,]$group = NA
snp_PCA_biplot_group = ggplot(data=snp_scores, aes(x=PC1, y=PC2)) +
  geom_hline(yintercept=0, colour="gray65") +
  geom_vline(xintercept=0, colour="gray65") +
  geom_mark_ellipse(data=subset(snp_scores, !is.na(group)), expand = 0.01, aes(color=group), alpha = 0, 
                    linetype = 'dashed', size = 0.4, show.legend = F) +
  geom_text(aes(label=rownames(snp_scores)),
            colour=NAM_colors, alpha=0.8, size=5) + 
  scale_color_manual(values=c('#53CACE', '#F5958F')) +
  xlab(paste("PC1 (", snp_PCvar[1], ")", sep="")) + 
  ylab(paste("PC2 (", snp_PCvar[2], ")", sep="")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14),
        legend.title=element_blank())
snp_PCA_biplot_group


###################################
#### BUSCO tree #####
###################################

library('ggtree')
busco_tree <- read.tree(file = 'busco-tree.nwk')
busco_tree$edge.length[is.na(busco_tree$edge.length)] = 0 #convert NAs to 0
busco_tree$tip.label = busco_tree$tip.label %>% sub("Oh7b", "Oh7B", .) %>% 
  sub("IL14H", "Il14H", .) %>% sub("MS71", "Ms71", .) %>% sub("Mo18*", "Mo18W", .)
tree_cols = c(NAM_colors, PARVIGLUMIS='black')
names(tree_cols) = c(NAM_colors, PARVIGLUMIS='black')
busco_tree = phytools::reroot(busco_tree, node.number=13) #use busco_tree$tip.label to find the node.number, read help
busco_tree = full_join(busco_tree, 
                       data.frame(label = names(c(NAM_colors, PARVIGLUMIS='black')), 
                                  stat = c(NAM_colors, PARVIGLUMIS='black')), by = 'label')
busco_tree_p <- ggtree(busco_tree, layout = 'rectangular', 
                       branch.length = 'none') + 
  geom_tiplab(aes(color = stat), size = 5, hjust = 0) + geom_rootpoint() +
  # geom_tippoint() +
  scale_color_manual(values = tree_cols) +
  theme(legend.position = "none") + xlim(NA, 15)
busco_tree_p


###################################
#### Most varying TE families #####
###################################
## calculate family variations
require(goeveg)
Var = apply(fam, 2, var, na.rm=TRUE)
class(Var)
str(Var)

## get top 10 varying TE families
rank = 10
var_rank = Var[order(-Var)][1:rank]
fam[,names(var_rank)]
topVar = data.frame(fam[names(var_rank)], Genome = rownames(fam))
str(topVar)
names(topVar)

# change the order of Genomes
new_order = c(NAM_ID[11:13], NAM_ID[1:10], NAM_ID[14:26])
topVar$Genome = factor(topVar$Genome, levels=new_order)
varlist = colnames(topVar)[1:10]

# gather plotting data
topVar_gather = topVar %>% gather(key=variable, value=value, names(topVar)[1:10])
topVar_gather$x=as.numeric(factor(topVar_gather$Genome)) #get the actual order of genome for plotting
str(topVar_gather)

# normalize the size
normalize <- function(x) {x / sqrt(sum(x^2))}
topVar_nor = data.frame(apply(topVar[,1:10], 2, scale), Genome = topVar$Genome) #standardization to mean=0, var=1
str(topVar_nor)

# gather plotting data
topVar_nor_gather = topVar_nor %>% gather(variable, value, names(topVar_nor)[1:10])
topVar_nor_gather$x = as.numeric(factor(topVar_nor_gather$Genome)) #get the actual order of genome for plotting
str(topVar_nor_gather)

# plot normalized family size with pvalue
boot.t.test(topVar_nor_gather[topVar_nor_gather$Genome %in% tropical, ]$value,
            topVar_nor_gather[topVar_nor_gather$Genome %in% non_tropical, ]$value)
2*pt(-abs(13.885),df=217.38) #actual p
topVar_nor_p = ggplot(topVar_nor_gather, aes(x = Genome, y = value)) + 
  geom_boxplot() + 
  geom_point(aes(color=variable), position = position_jitter(width = 0.15), size = 1) +
  geom_signif(comparisons=list(c("HP301", "Oh7B")), annotations="", y_position = 1.5, 
              tip_length = 0, color = '#53CACE', size = 0.5) +
  geom_signif(comparisons=list(c("CML52", "Tzi8")), annotations="", y_position = 2.5, 
              tip_length = 0, color = '#F5958F', size = 0.5) +
  geom_signif(annotations="P = 8.4e-32", y_position = 3.1, xmin=5.5, xmax=20, size = 0.5, 
              vjust = -0.5, tip_length = c(0.21, 0.03)) + 
  labs(x =" ", y = "Standardized family size",size=14, face="bold") +
  ylim(-4.4, 3.5) +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(color = NAM_colors[new_order], size=10, angle=35, vjust=1, hjust=0.95),
        legend.text = element_text(size=8),
        legend.title=element_blank()) +
  theme(legend.position = c(0.5, 0.1), legend.direction = "horizontal", 
        legend.text = element_text(size=5.5),
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.3,"cm")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
topVar_nor_p


###################################
#### Top 50 largest TE families #####
###################################

# Plot out largest TE families and variation in NAMs
str(fam)
rank = 50
rownames(fam) #genome name
std <- function(x) sd(x)/sqrt(length(x)) #standard error
fam_size = data.frame(TE = colnames(fam), 
                      sd = apply(fam, 2, sd, na.rm=TRUE), 
                      se = apply(fam, 2, std),
                      mean = apply(fam, 2, mean, na.rm=TRUE),
                      Min = apply(fam, 2, min, na.rm=TRUE), 
                      Max = apply(fam, 2, max, na.rm=TRUE))
str(fam_size)
topSize = fam_size[order(-fam_size$mean),][1:rank,]
topSize$TE = factor(topSize$TE, levels = topSize$TE) #control the plotting order
head(topSize)
topSize_100 = fam_size[order(-fam_size$mean),][1:100,]
sum(topSize$mean) #top 50
sum(topSize_100$mean) # top 100
sum(fam_size$mean) #size of all TEs

# top size plot, error bar SD
topSize_plot_SD = ggplot(topSize, aes(x=TE, y=mean)) + 
  geom_bar(stat="identity", fill = rgb(27,168,241, max=255)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  labs(x ="TE family", y = "Family size (Mb)", size=12, face="bold")+
  theme(axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=8, angle=50,vjust=1,hjust=0.95),
        plot.margin = margin(0.2, 0.2, 0.2, 0.7, "cm")
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
topSize_plot_SD

##############################
#### Compile intact LTRs  ####
##############################
#read 26 genome LTR data and calculate LTR age
LTR_len = read.table('NAM.26.intact.LTR.len.info', header = F)
names(LTR_len) = c('genome', 'TE_id', 'id', 'len', 'ltrlen', 'intlen')
LTR = read.table('NAM.26.intact.LTR.list', header = F, fill = T)
notLTR = read.table('NAM.26.intact.not_LTR.list', header = F, fill = T)
str(LTR)
colnames(LTR) = c('genome', 'chr', 'supfam', 'classification', 'from', 'to', 'strand', 'id', 'SO', 'motif', 'tsd', 'iden')
colnames(notLTR) = c('genome', 'chr', 'supfam', 'classification', 'from', 'to', 'strand', 'id', 'SO')
LTR$age = (1-LTR$iden)/(2*3.3e-8)/1000
LTR$len = LTR$to - LTR$from + 1
notLTR$len = notLTR$to - notLTR$from + 1
LTR$TE_id = paste(LTR$chr, ':', LTR$from, '..', LTR$to, sep = '')
LTR$ltrlen = LTR_len$ltrlen[match(LTR$TE_id, LTR_len$TE_id)]
LTR$intlen = LTR_len$intlen[match(LTR$TE_id, LTR_len$TE_id)]
LTR$group = LTR$genome
LTR[LTR$group %in% tropical,]$group = 'trop'
LTR[LTR$group %in% non_tropical,]$group = 'temp'
LTR[!LTR$group %in% c('temp', 'trop'),]$group = 'admx'
rm(LTR_len)
head(LTR)
head(notLTR)

# fix genome names
LTR$genome = sub("Oh7b", "Oh7B", LTR$genome)
LTR$id = sub("Oh7b", "Oh7B", LTR$id)
notLTR$genome = sub("Oh7b", "Oh7B", notLTR$genome)
notLTR$id = sub("Oh7b", "Oh7B", notLTR$id)
unique(LTR$genome)
unique(notLTR$genome)

# count
length(LTR$id)
sum(str_detect(LTR$id, ":"))
sum(str_detect(LTR$id, ":", negate = T))
table(LTR$genome)
mean(table(LTR$genome)) #53678.5
sd(table(LTR$genome)) #1347.759

# a function to modify facet labels
label_facet <- function(original_var, custom_name){
  lev <- levels(as.factor(original_var))
  custom_name = round(custom_name, 2)
  lab <- paste0(lev, " (", custom_name, " Mb)")
  names(lab) <- lev
  return(lab)  
}

LTR_mean_age = aggregate(sizelist_LTR$age, list(sizelist_LTR$id), mean)
LTR_mean_age$status_class = Fam_count$status_class[match(LTR_mean_age$Group.1, Fam_count$TE)]
LTR_mean_age = LTR_mean_age[with(LTR_mean_age, order(status_class, x)),]


##############################
#### Annotate intact LTRs ####
##############################
LTR_sort = read.table('NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.intact.LTR.sort', header = F, sep = '\t')
names(LTR_sort) = c('genome', 'TE_id', 'id', 'Order', 'Superfamily', 'Clade', 'Complete', 'Strand', 'Domains', 'domain_ct')
str(LTR_sort)

# remove the B73_ab10 and fix genome names
LTR_sort = LTR_sort[!LTR_sort$genome %in% c("B73_AB10", "B73Ab10"),]
LTR_sort$genome = sub("Oh7b", "Oh7B", LTR_sort$genome)
unique(LTR_sort$genome)

# plot domain count freq
LTR_sort[LTR_sort$Complete %in% c('none', 'unknown'),]$Complete = 'no'
LTR_sort_domct_plot = ggplot(LTR_sort, aes(domain_ct)) + geom_histogram() + 
  facet_wrap(~ Complete) + #title('Domain completeness') +
  labs(x='Domain numbers', y='Element count') + theme_bw()
LTR_sort_domct_plot

##############################
#### Trop/Temp LTRs variation ####
##############################
# find LTR fams with size sig. diff in trop and temp
str(fam)
fam_nontrop = fam[rownames(fam) %in% NAM_nontrop,]
fam_trop = fam[rownames(fam) %in% NAM_trop,]
test_trop_temp <- vector("list", 1)
for (j in seq(ncol(fam))){
  temp = fam_nontrop[,j]
  trop = fam_trop[,j]
  if (sd(temp-trop) == 0){
    temp = jitter(temp)
  }
  test_trop_temp[[j]] = t.test(temp, trop)
}
test_trop_temp[[1]]
test_trop_temp = data.frame(matrix(unlist(test_trop_temp), nrow=length(test_trop_temp), byrow=T))
test_trop_temp[,1:9] = lapply(test_trop_temp[,1:9], function(x) {as.numeric(x)})
names(test_trop_temp)[names(test_trop_temp) == 'X3'] = 'pval'
names(test_trop_temp)[names(test_trop_temp) == 'X6'] = 'temp_mean'
names(test_trop_temp)[names(test_trop_temp) == 'X7'] = 'trop_mean'
test_trop_temp$diff = test_trop_temp$trop_mean - test_trop_temp$temp_mean
test_trop_temp$mean = (test_trop_temp$trop_mean*13 + test_trop_temp$temp_mean*10)/23
rownames(test_trop_temp) = sub('_INT|_LTR', '', colnames(fam))
test_trop_temp$id = rownames(test_trop_temp)
str(test_trop_temp)

# stats
dim(subset(test_trop_temp, temp_mean==0))
dim(subset(test_trop_temp, trop_mean==0))
sum(subset(test_trop_temp, temp_mean==0)$diff)
sum(subset(test_trop_temp, trop_mean==0)$diff)

# find LTR fams with intact size sig. diff b/t trop and temp
str(fam_int)
rowSums(fam_int[,colnames(fam_int) %in% Fam_count$TE])
sum(str_detect(names(fam_int), ':', negate = T), na.rm = T)
fam_int_pan = fam_int[, names(fam_int) %in% Fam_count$TE]
fam_int_nontrop = fam_int_pan[rownames(fam_int_pan) %in% NAM_nontrop,]
fam_int_trop = fam_int_pan[rownames(fam_int_pan) %in% NAM_trop,]
test_trop_temp_int <- vector("list", 1)
for (j in seq(ncol(fam_int_pan))){
  temp = fam_int_nontrop[,j]
  trop = fam_int_trop[,j]
  if (sd(temp-trop) == 0){
    temp = jitter(temp)
  }
  test_trop_temp_int[[j]] = t.test(temp, trop)
}
test_trop_temp_int[[1]]
test_trop_temp_int = data.frame(matrix(unlist(test_trop_temp_int), nrow=length(test_trop_temp_int), byrow=T))
test_trop_temp_int[,1:9] = lapply(test_trop_temp_int[,1:9], function(x) {as.numeric(x)})
names(test_trop_temp_int)[names(test_trop_temp_int) == 'X3'] = 'pval'
names(test_trop_temp_int)[names(test_trop_temp_int) == 'X6'] = 'temp_mean'
names(test_trop_temp_int)[names(test_trop_temp_int) == 'X7'] = 'trop_mean'
test_trop_temp_int$diff = test_trop_temp_int$trop_mean - test_trop_temp_int$temp_mean
test_trop_temp_int$mean = (test_trop_temp_int$trop_mean*13 + test_trop_temp_int$temp_mean*10)/23
rownames(test_trop_temp_int) = sub('_INT|_LTR', '', colnames(fam_int_pan))
test_trop_temp_int$id = rownames(test_trop_temp_int)
test_trop_temp_int$category = Fam_count$category[match(test_trop_temp_int$id, Fam_count$TE)]
str(test_trop_temp_int)

# stats
sum(test_trop_temp_int$diff)
sum(test_trop_temp_int$diff)/sum(test_trop_temp$diff) #46.3%
sum(fam_int)/26
554.636/1796.746 #0.3086892
fisher.test(rbind(c(19.60124,22.76089), c(554.636, 1242.11))) #p-value = 0.04408
tapply(test_trop_temp_int$diff, test_trop_temp_int$category, sum)
dim(subset(test_trop_temp_int, pval<0.05))
table(subset(test_trop_temp_int, pval<0.05)$category)
Fam_count$int_diff_pval = test_trop_temp_int$pval[match(Fam_count$TE, test_trop_temp_int$id)]
Fam_count$trop_temp_diff_int = test_trop_temp_int$diff[match(Fam_count$TE, test_trop_temp_int$id)]
test_trop_temp_int$class = class$Supfam[match(test_trop_temp_int$id, sub('_INT|_LTR', '', class$TE))]
test_trop_temp_int = test_trop_temp_int[!is.na(test_trop_temp_int$class),]
table(test_trop_temp_int$class)
rowSums(fam_int)

# stat LTR fams (more stats in the solo-intact section)
test_trop_temp[test_trop_temp$diff>1,]
trop_temp_difflist = colnames(fam)[test_trop_temp$diff>1]
fam[colnames(fam)[test_trop_temp$diff>1]]
sum(test_trop_temp$diff)
sum(test_trop_temp[test_trop_temp$diff<0,]$diff)
sum(test_trop_temp[test_trop_temp$diff>0,]$diff)

# make groups
cutoff = 0.025
test_trop_temp = test_trop_temp %>%
  mutate(diff_group = case_when(diff < -cutoff ~ paste('(-inf, ', -cutoff, ']', sep = ''),
                                diff >= -cutoff & diff <= cutoff ~ paste('(', -cutoff, ', ', cutoff, ')', sep = ''),
                                diff > cutoff ~ paste('[', cutoff, ', inf)', sep = '')))
test_trop_temp$diff_group = as.factor(test_trop_temp$diff_group)
test_trop_temp$diff_group <- factor(test_trop_temp$diff_group, 
                                    levels = c(paste('(-inf, ', -cutoff, ']', sep = ''),
                                               paste('(', -cutoff, ', ', cutoff, ')', sep = ''),
                                               paste('[', cutoff, ', inf)', sep = '')))
trop_temp_diff_rank_plot4 = ggplot() + 
  geom_violin(data = test_trop_temp, aes(diff_group, diff), fill = 'grey', scale = "width") +
  geom_text(aes(label = paste('n = ', format(dim(subset(test_trop_temp, diff < -cutoff))[1], big.mark=","), sep = '')), x = 1, y = -0.6, size = 5) +
  geom_text(aes(label = paste('n = ', format(dim(subset(test_trop_temp, diff >= -cutoff & diff <= cutoff))[1], big.mark=","), sep = '')), x = 2, y = -0.6, size = 5) +
  geom_text(aes(label = paste('n = ', format(dim(subset(test_trop_temp, diff > cutoff))[1], big.mark=","), sep = '')), x = 3, y = -0.6, size = 5) +
  labs(y = 'Trop - Temp family size (Mb)', x = 'Temp. larger               No diff.               Trop. larger') +
  ylim(-0.7, 4.7) +
  theme(axis.title.x = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        axis.text.x = element_text(size=10, vjust = -3),
        axis.text.y = element_text(size=10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
trop_temp_diff_rank_plot4

# make groups for fam size
cutoff = 0.1
fam_size = fam_size %>%
  mutate(size_group = case_when(mean <= cutoff ~ paste('[0, ', cutoff, ']', sep = ''),
                                mean > cutoff ~ paste('(', cutoff, ', inf)', sep = '')))
fam_size$size_group = as.factor(fam_size$size_group)
fam_size$size_group <- factor(fam_size$size_group, 
                              levels = c(paste('[0, ', cutoff, ']', sep = ''),
                                         paste('(', cutoff, ', inf)', sep = '')))
fam_size_rank_plot = ggplot() + 
  geom_violin(data = fam_size, aes(size_group, mean), scale = "width") +
  geom_text(aes(label = paste('n = ', 27228 - dim(subset(fam_size, mean > cutoff))[1], sep = '')), x = 1.25, y = 10, size = 2.8) +
  geom_text(aes(label = paste('n = ', dim(subset(fam_size, mean > cutoff))[1], sep = '')), x = 2.25, y = 10, size = 2.8) +
  labs(y = 'Family size (Mb)', x = '') +
  theme(axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10, vjust = -1),
        axis.text.y = element_text(size=10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
fam_size_rank_plot

test_trop_temp$class = sub('LINE.*', 'LINE', test_trop_temp$class)
trop_temp_size_cmb_plot = ggplot(test_trop_temp, aes(x = class, y = diff, fill = class, color = class)) +
  geom_bar(stat = 'identity', size = 0.05) +
  geom_hline(yintercept=0, size = 0.3, color = 'black') +
  scale_fill_manual(values=TE_colors3) +
  scale_color_manual(values=TE_colors3) +
  ylab("Trop - Temp TE size (Mb)") +
  scale_y_continuous(breaks = c(-5, 0, 10, 20, 30), limits = c(-6, 31)) +
  theme(legend.position="none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=12, angle=35, vjust=1, hjust=0.95)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
trop_temp_size_cmb_plot

test_trop_temp_int$class = sub('LINE.*', 'LINE', test_trop_temp_int$class)
trop_temp_size_int_cmb_plot = ggplot(test_trop_temp_int, aes(x = class, y = diff, fill = class, color = class)) +
  geom_bar(stat = 'identity', size = 0.05) +
  geom_hline(yintercept=0, size = 0.3, color = 'black') +
  scale_fill_manual(values=TE_colors3) +
  scale_color_manual(values=TE_colors3) +
  ylab("Trop - Temp Intact size (Mb)") +
  scale_y_continuous(breaks = c(-5, 0, 5, 10, 15), limits = c(-6, 15)) +
  theme(legend.position="none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=12, angle=35, vjust=1, hjust=0.95)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
trop_temp_size_int_cmb_plot

# size summary
sum(test_trop_temp$trop_mean) - sum(test_trop_temp$temp_mean)
sum(test_trop_temp[str_detect(test_trop_temp$class, 'LTR'),]$diff)
sum(test_trop_temp[str_detect(test_trop_temp$class, 'TIR'),]$diff)
sum(test_trop_temp[str_detect(test_trop_temp$class, 'CACTA'),]$diff)
sum(test_trop_temp[str_detect(test_trop_temp$class, 'Helitron'),]$diff)
sum(test_trop_temp[str_detect(test_trop_temp$class, 'LINE'),]$diff)


###################################
#### Classify each LTR family  ####
###################################
str(LTR) #all intact LTRs in NAM

# draw the classification scheme
act2 = ggplot(data = data.frame(x = c(-3000, 3000)), aes(x)) +
  stat_function(fun = dnorm, n = 1001, args = list(mean = 0, sd = 500)) + labs(y="", x='') +
  scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(200, 3000)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab('Y(max) = Y(0)') + ylab('Active') +
  theme(axis.title.x = element_text(vjust = -1), axis.title.y = element_text(size = 12, face = 'bold'))
phd3 = ggplot(data = data.frame(x = c(-3000, 3000)), aes(x)) +
  geom_hline(yintercept = 0.0001, linetype = 'dashed', color = 'red') +
  stat_function(fun = dnorm, n = 1001, args = list(mean = 0, sd = 500)) + labs(y="", x='') +
  scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(-300, 3000)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab('Y(max) > Y(0); Y(0) >= 5% total') + ylab('Phasing down') +
  theme(axis.title.x = element_text(vjust = -1), axis.title.y = element_text(size = 12, face = 'bold'))
ina2 = ggplot(data = data.frame(x = c(-3000, 3000)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 500)) + labs(y="", x='') +
  scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(-1500, 3000)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab('Y(max) > Y(0); Y(0) < 5% total') +
  theme(axis.title.x = element_text(vjust = -1, hjust = -1),
        plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 1.1, unit = "in"))


## a function to classify LTR fam type
classify_fam = function(x) {
  modes = NULL #distribution mode
  type = NULL #family type
  
  # get data and basic stats
  plot = ggplot(x, aes(1-iden)) + geom_histogram(binwidth=0.002, center=0.001) #type is sensitive to binwidth
  y = ggplot_build(plot)$data[[1]]$count
  xmin = ggplot_build(plot)$data[[1]]$xmin[1]
  if (xmin == 0){
    density0 = y[1]/sum(y)
    maxy = which.max(y)
  } else {
    density0 = 0
    maxy = 'NA'
  }
  
  stats = c(density0, maxy, sum(y))
  names(stats) = c('density0', 'max_y', 'sum_y')
  
  # determine LTR activities
  if (maxy == 1 & xmin == 0){
    type = "active"
  } else if (density0 >= 0.05 & xmin == 0){
    type = "phasing down"
  } else {
    type = "inactive"
  }
  
  # determine the number of peaks
  if (length(y) < 5){
    modes = c(modes,maxy)
  } else {
    if (y[1] >= y[2] & y[1] >= y[3] & y[1]/sum(y) >=0.05) {
      modes = c(modes,1)
    }
    if (y[2] >= y[1] & (y[2] >= y[3] & y[2] >= y[4]) & y[2]/sum(y) >=0.05) {
      modes = c(modes,2)
    }
    for ( i in 3:(length(y)-2) ){
      if ( (y[i] >= y[i-2] & y[i] >= y[i-1]) & (y[i] >= y[i+1] & y[i] >= y[i+2]) & y[i]/sum(y) >=0.05 ) {
        modes <- c(modes,i)
      }
    }
  }
  if ( length(modes) == 0 ) {
    modes = 'This is a monotonic distribution'
  }
  
  results = list(type, modes, stats)
  return(results)
}

## classify shared LTR families
str(LTR_shared)
LTR_fam_list = unique(LTR_shared$id) #use LTR families shared >= 2 genomes, total 96815 fams
LTR_fam_list = LTR_fam_list[!str_detect(LTR_fam_list, ":")] #remove genome specific LTRs
LTR_fam_info = read.csv(text="id,class,modal,density0,max_y,sum_y") #create an empty df
for ( i in 1:length(LTR_fam_list) ){
  print(LTR_fam_list[i])
  this_fam = LTR[LTR$id==LTR_fam_list[i],]
  if (length(this_fam$age)>1){
    info = classify_fam(this_fam)
    fam_info = data.frame(id = LTR_fam_list[i],
                          class = info[[1]], 
                          modal = length(info[[2]]), 
                          density0 = info[[3]][1],
                          max_y = info[[3]][2],
                          sum_y = info[[3]][3]
    )
    LTR_fam_info = rbind(LTR_fam_info, fam_info)
  }
}
str(LTR_fam_info)

# summarize
str(LTR_shared)
LTR_fam_info$sum_y = as.numeric(LTR_fam_info$sum_y)
LTR_fam_info$len_mean = LTR_shared$size_sum[match(LTR_fam_info$id, LTR_shared$id)]/LTR_fam_info$sum_y #mean LTR len of each fam in NAM
LTR_fam_info$age_mean = LTR_shared$age_sum[match(LTR_fam_info$id, LTR_shared$id)]/LTR_fam_info$sum_y #mean LTR age of each fam in NAM
LTR_fam_info$size_mean = LTR_shared$size_sum[match(LTR_fam_info$id, LTR_shared$id)]/26/1000000 #mean LTR fam size in each NAM
table(LTR_fam_info$modal)
table(LTR_fam_info$class)

# mean family count in each NAM
mean(LTR_fam_info[LTR_fam_info$class=='active',]$sum_y)/26
mean(LTR_fam_info[LTR_fam_info$class=='phasing down',]$sum_y)/26
mean(LTR_fam_info[LTR_fam_info$class=='inactive',]$sum_y)/26

# mean family size in each NAM (Mb)
mean(LTR_fam_info[LTR_fam_info$class=='active',]$size_mean)
mean(LTR_fam_info[LTR_fam_info$class=='phasing down',]$size_mean)
mean(LTR_fam_info[LTR_fam_info$class=='inactive',]$size_mean)

## a simplified function to classify LTR fam type
classify_fam_simp = function(x) {
  modes = NULL #distribution mode
  type = NULL #family type
  
  # get data and basic stats
  plot = ggplot(x, aes(1-iden)) + geom_histogram(binwidth=0.002, center=0.001) #type is sensitive to binwidth
  y = ggplot_build(plot)$data[[1]]$count
  xmin = ggplot_build(plot)$data[[1]]$xmin[1]
  if (xmin == 0){
    density0 = y[1]/sum(y)
    maxy = which.max(y)
  } else {
    density0 = 0
    maxy = 'NA'
  }
  
  # determine LTR activities
  if (maxy == 1 & xmin == 0){
    type = "active" # == 'Young'
  } else if (density0 >= 0.05 & xmin == 0){
    type = "phasing down" # == 'Moderate'
  } else {
    type = "inactive" # == 'Old'
  }
  return(type)
}


#########################
#### solo/intact LTR ####
#########################

# read solo LTR info
solo = read.table('NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.solo', header = F)
colnames(solo) = c('genome', 'chr', 'from', 'to', 'solo_id', 'id', 'cov')
solo$id = sub('_LTR', '', solo$id) #remove _LTR string in id
solo = solo[solo$genome != 'B73Ab10' & solo$genome != 'B73_AB10', ] #rm Ab10
solo$genome = as.character(solo$genome) #rm Ab10 factors
str(solo)

# fix genome names
solo$genome = sub("Oh7b", "Oh7B", solo$genome)
solo$id = sub("Oh7b", "Oh7B", solo$id)
solo = solo[!str_detect(solo$genome, 'alt|scf'),]
unique(solo$genome)

# count
head(sort(table(solo$id), decreasing = T))
head(sort(table(LTR$id), decreasing = T))
table(solo$id)/table(LTR$id)
length(table(solo$id))
length(unique(LTR[!LTR$id %in% str_subset(LTR$id, ":"),]$id))
length(unique(solo[!solo$id %in% str_subset(solo$id, ":"),]$id))
str_subset(unique(LTR[!LTR$id %in% str_subset(LTR$id, ":"),]$id), 'TE_00019874')
str_subset(unique(solo[!solo$id %in% str_subset(solo$id, ":"),]$id), 'TE_00018037')

# calculate solo/intact ratio for each fam in each NAM
LTRlist = unique(LTR[!LTR$id %in% str_subset(LTR$id, ":"),]$id)
solo_intact = data.frame(id = LTRlist)
for (i in 1:length(NAM_ID)){
  print (NAM_ID[i])
  intact_i = LTR[LTR$genome == NAM_ID[i] & LTR$id %in% LTRlist,]
  solo_i = solo[solo$genome == NAM_ID[i] & solo$id %in% LTRlist,]
  intact_count = table(as.character(intact_i$id))
  solo_count = table(as.character(solo_i$id))
  solo_intact_i = c()
  for (j in 1:length(LTRlist)){
    intact_cj = intact_count[LTRlist[j]]
    solo_cj = solo_count[LTRlist[j]]
    if (is.na(intact_cj)){
      intact_cj = 1
    }
    if (is.na(solo_cj)){
      solo_cj = 0
    }
    si = solo_cj/intact_cj
    solo_intact_i = c(solo_intact_i, si)
  }
  solo_intact[ , ncol(solo_intact) + 1] = solo_intact_i
  colnames(solo_intact)[ncol(solo_intact)] = NAM_ID[i]
  # solo_intact[[i]] = solo_intact_i
}
str(solo_intact)


####################################
#### Process UM CHG LTRs  ####
####################################
## get UMR originating from lLTR
#UMR_sum = read.table('NAM26.intact.LTR.UMR.sum', header = F)
UMR_sum = read.table('./methylation/NAM26.intact.lLTR.UMR.nostart', header = F)
UMR_sum$V9 = NULL
UMR_sum$V13 = NULL
UMR_sum$V14 = NULL
names(UMR_sum) = c('genome', 'geneloca', 'chr', 'from', 'to', 'TE_id', 'strand', 'iden', 
                   'umfrom', 'umto', 'umid', 'umlen')
str(UMR_sum)
UMR_sum$genome = UMR_sum$genome %>% sub("OH7b", "Oh7B", ignore.case = T, .) %>% 
  sub("IL14H", "Il14H", ignore.case = T, .) %>% 
  sub("MS71", "Ms71", ignore.case = T, .) %>% 
  sub("OH43", "Oh43", ignore.case = T, .) %>%
  sub("TZi8", "Tzi8", ignore.case = T, .) %>%
  sub("MS37W", "M37W", ignore.case = T, .) %>%
  sub("TX303", "Tx303", ignore.case = T, .)
UMR_sum = UMR_sum %>% mutate(UM_dist = ifelse(strand == '+', umfrom - from, to - umto))
t.test(subset(UMR_sum, genome %in% NAM_trop)$iden, 
       subset(UMR_sum, genome %in% NAM_nontrop)$iden)
t.test(subset(UMR_sum, genome %in% NAM_trop)$umlen, 
       subset(UMR_sum, genome %in% NAM_nontrop)$umlen)
t.test(UMR_sum$umlen ~ UMR_sum$geneloca)

LTR$umlen = UMR_sum$umlen[match(LTR$TE_id, sub('.*_', '', UMR_sum$TE_id))]
LTR = LTR %>% mutate(CHG_UM = case_when(umlen >= 0 ~ 'yes', TRUE ~ 'no'))
LTR$UM_dist = UMR_sum$UM_dist[match(LTR$TE_id, sub('.*_', '', UMR_sum$TE_id))]
tapply(UMR_sum$umlen, UMR_sum$strand, mean)
tapply(subset(LTR, umlen>0 & geneloca != 'NA' & classification == 'LTR/Copia')$umlen,
       subset(LTR, umlen>0 & geneloca != 'NA' & classification == 'LTR/Copia')$geneloca, mean)
tapply(subset(LTR, umlen>0 & geneloca != 'NA' & classification == 'LTR/Gypsy')$umlen,
       subset(LTR, umlen>0 & geneloca != 'NA' & classification == 'LTR/Gypsy')$geneloca, mean)
cor.test(UMR_sum$umlen, UMR_sum$iden) #p-value = 0.08494
t.test(UMR_sum$umlen ~ UMR_sum$geneloca) #p-value = 0.03396, inside longer
t.test(subset(LTR, classification != 'LTR/unknown')$umlen ~ #p-value < 2.2e-16, Gypsy longer
         subset(LTR, classification != 'LTR/unknown')$classification)
dim(subset(LTR, category == 'Tropical amplification' & umlen >0))[1]/dim(subset(LTR, umlen >0))[1]
test = as.data.frame(sort(table(subset(LTR, umlen > 0 & str_detect(id, ':', negate = T))$id), decreasing = T))


# get certain LTRs
# xilon_diguus_AC203313_7774, trop amp with biggest diff
UM_LTR_topA_list = paste(subset(LTR, umlen > 0 & id == 'xilon_diguus_AC203313_7774' & geneloca == 'outside')$genome, 
                         subset(LTR, umlen > 0 & id == 'xilon_diguus_AC203313_7774' & geneloca == 'outside')$TE_id, sep = '_')
UM_LTR_topA_p1 = ggplot(subset(LTR, umlen > 0 & id == 'xilon_diguus_AC203313_7774' & geneloca == 'outside')) +
  geom_histogram(aes(ltrlen)) + labs(y = 'Count', x = 'LTR length') + theme_bw()
UM_LTR_topA_p2 = ggplot(subset(LTR, umlen > 0 & id == 'xilon_diguus_AC203313_7774' & geneloca == 'outside')) +
  geom_histogram(aes(umlen)) + labs(y = 'Count', x = 'UMR length') + theme_bw()
UM_LTR_topA_p3 = ggplot(subset(UMR_sum, sub('.*_', '', TE_id) %in% 
                                 subset(LTR, umlen > 0 & id == 'xilon_diguus_AC203313_7774' & geneloca == 'outside')$TE_id)) +
  geom_histogram(aes(UM_dist)) + labs(y = 'Count', x = 'UMR distance') + theme_bw()
summary(subset(UMR_sum, sub('.*_', '', TE_id) %in% 
                 subset(LTR, umlen > 0 & id == 'xilon_diguus_AC203313_7774' & geneloca == 'outside')$TE_id)$lltr_dist)
UM_LTR_topA_p1 | UM_LTR_topA_p2 | UM_LTR_topA_p3
write.table(UM_LTR_topA_list, file = "UM_LTR_topA_list.txt", append = FALSE, quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)

# ji_AC215728_13156, temp rmv with biggest diff
UM_LTR_topB_list = subset(LTR, umlen > 0 & id == 'ji_AC215728_13156' & geneloca == 'outside')$TE_id
UM_LTR_topB_p1 = ggplot(subset(LTR, umlen > 0 & id == 'ji_AC215728_13156' & geneloca == 'outside')) +
  geom_histogram(aes(ltrlen)) + labs(y = 'Count', x = 'LTR length') + theme_bw()
UM_LTR_topB_p2 = ggplot(subset(LTR, umlen > 0 & id == 'ji_AC215728_13156' & geneloca == 'outside')) +
  geom_histogram(aes(umlen)) + labs(y = 'Count', x = 'UMR length') + theme_bw()
UM_LTR_topB_p3 = ggplot(subset(UMR_sum, sub('.*_', '', TE_id) %in% 
                                 subset(LTR, umlen > 0 & id == 'ji_AC215728_13156' & geneloca == 'outside')$TE_id)) +
  geom_histogram(aes(UM_dist)) + labs(y = 'Count', x = 'UMR distance') + theme_bw()
summary(subset(UMR_sum, sub('.*_', '', TE_id) %in% 
                 subset(LTR, umlen > 0 & id == 'ji_AC215728_13156' & geneloca == 'outside')$TE_id)$lltr_dist)
UM_LTR_topB_p1 + UM_LTR_topB_p2 + UM_LTR_topB_p3 + plot_annotation(title = 'ji_AC215728_13156')

# uwum_AC213069_12092, trop amp, most umr
UM_LTR_topC_list = subset(LTR, umlen > 0 & id == 'uwum_AC213069_12092' & geneloca == 'outside')$TE_id
UM_LTR_topC_p1 = ggplot(subset(LTR, umlen > 0 & id == 'uwum_AC213069_12092' & geneloca == 'outside')) +
  geom_histogram(aes(ltrlen)) + labs(y = 'Count', x = 'LTR length') + theme_bw()
UM_LTR_topC_p2 = ggplot(subset(LTR, umlen > 0 & id == 'uwum_AC213069_12092' & geneloca == 'outside')) +
  geom_histogram(aes(umlen)) + labs(y = 'Count', x = 'UMR length') + theme_bw()
UM_LTR_topC_p3 = ggplot(subset(UMR_sum, sub('.*_', '', TE_id) %in% 
                                 subset(LTR, umlen > 0 & id == 'uwum_AC213069_12092' & geneloca == 'outside')$TE_id)) +
  geom_histogram(aes(UM_dist)) + labs(y = 'Count', x = 'UMR distance') + theme_bw()
summary(subset(UMR_sum, sub('.*_', '', TE_id) %in% 
                 subset(LTR, umlen > 0 & id == 'uwum_AC213069_12092' & geneloca == 'outside')$TE_id)$lltr_dist)
UM_LTR_topC_p1 | UM_LTR_topC_p2 | UM_LTR_topC_p3

# CRM2_7577nt, trop amp, known active
UM_LTR_topD_list = paste(subset(LTR, umlen > 0 & id == 'CRM2_7577nt' & geneloca == 'outside' & strand == '+')$genome, 
                         subset(LTR, umlen > 0 & id == 'CRM2_7577nt' & geneloca == 'outside' & strand == '+')$TE_id, sep = '_')
UM_LTR_topD_p1 = ggplot(subset(LTR, umlen > 0 & id == 'CRM2_7577nt' & geneloca == 'outside')) +
  geom_histogram(aes(ltrlen)) + labs(y = 'Count', x = 'LTR length') + theme_bw()
UM_LTR_topD_p2 = ggplot(subset(LTR, umlen > 0 & id == 'CRM2_7577nt' & geneloca == 'outside')) +
  geom_histogram(aes(umlen)) + labs(y = 'Count', x = 'UMR length') + theme_bw()
UM_LTR_topD_p3 = ggplot(subset(UMR_sum, sub('.*_', '', TE_id) %in% 
                                 subset(LTR, umlen > 0 & id == 'CRM2_7577nt' & geneloca == 'outside')$TE_id)) +
  geom_histogram(aes(UM_dist)) + labs(y = 'Count', x = 'UMR distance') + theme_bw()
summary(subset(UMR_sum, sub('.*_', '', TE_id) %in% 
                 subset(LTR, umlen > 0 & id == 'CRM2_7577nt' & geneloca == 'outside')$TE_id)$lltr_dist)
write.table(UM_LTR_topD_list, file = "UM_LTR_topD_list.txt", append = FALSE, quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)


# plot
LTR_umlen_plot = ggplot(subset(LTR, umlen>0 & geneloca != 'NA' & classification != 'LTR/unknown')) + 
  geom_density(aes(umlen, color = geneloca)) +
  scale_color_manual(name='Position to genes', values = c("#0072B2", "#D55E00"), 
                     labels = c("Overlapping", "Outside")) +
  scale_x_continuous(trans='log10') + facet_wrap(~classification) + theme_bw() +
  labs(x='UM track length (bp)', y='Density') + theme(legend.position = c(0.9, 0.75))
LTR_umlen_plot

# aggregate
# count and add info to these intact LTRs
LTR_um_group = subset(LTR, umlen>0 & geneloca != 'NA') %>% group_by(group, id) %>% 
  summarize(count=n(), len = sum(to - from + 1), 
            id_group = unique(paste(id, group, sep = '_')),
            Supfam = names(sort(table(classification), decreasing = T)[1])) #get supfam from major consensus
LTR_um_group$supfam = LTR_genes$supfam[match(LTR_um_group$id, LTR_genes$id)]
LTR_um_group$Supfam = factor(LTR_um_group$Supfam, levels = c('LTR/Gypsy', 'LTR/Copia', 'LTR/unknown'))
LTR_um_group$condition = LTR_genes$condition[match(LTR_um_group$id, LTR_genes$id)]
LTR_um_group$status_class = LTR_genes$status_class[match(LTR_um_group$id, LTR_genes$id)] #trop temp together classified
LTR_um_group[is.na(LTR_um_group$status_class),]$status_class = 'unknown'
LTR_um_group$status_class = factor(LTR_um_group$status_class, levels = c('active', 'phasing down', 'inactive', 'unknown'))
LTR_um_group$group = factor(LTR_um_group$group, levels=c('trop', 'temp', 'admx'))
LTR_um_group$category = Fam_count$category[match(LTR_um_group$id, Fam_count$TE)]
LTR_um_group$trop_temp_diff = Fam_count$trop_temp_diff[match(LTR_um_group$id, Fam_count$TE)] 

# plot trop temp UM counts
LTR_um_group_count_plot2 = ggplot() + 
  geom_bar(data = subset(LTR_um_group, Supfam != 'LTR/unknown' & group == 'trop' & trop_temp_diff > 0 &
                           status_class != 'unknown'), stat = 'identity', alpha = 1, width = 0.35, 
           aes(as.numeric(status_class)+0.05, count, fill = category, color = '#1fbfc3', size = 2))  +
  geom_bar(data = subset(LTR_um_group, Supfam != 'LTR/unknown' & group == 'temp' & trop_temp_diff > 0 &
                           status_class != 'unknown'), stat = 'identity', alpha = 1, width = 0.35,
           aes(as.numeric(status_class)+0.45, count, fill = category, color = '#f57770', size = 2)) +
  geom_text(data = subset(LTR_um_group, group == 'trop' & status_class != 'unknown' & Supfam != 'LTR/unknown'), stat = "count", 
            aes(as.numeric(status_class)+0.05, count, label = ..count..), size = 2, y=0, vjust = 1.2) +
  geom_text(data = subset(LTR_um_group, group == 'temp' & status_class != 'unknown' & Supfam != 'LTR/unknown' ), stat = "count", 
            aes(as.numeric(status_class)+0.45, count, label = ..count..), size = 2, y=0, vjust = 1.2) +
  labs(x ="") + ylab("Unmethylated intact LTR-RTs (#)") +
  scale_size(range = c(0.2, 0.4), guide='none') +
  facet_grid(.~Supfam) +
  guides(color = guide_legend(override.aes = list(fill = "white"))) +
  scale_fill_manual(breaks = c('Temperate removal', 'Balanced (trop > temp)', 
                               'Tropical amplification', 'Drifting (trop > temp)'),
                    values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                    name = 'Classifications') +
  scale_color_discrete(name = '', breaks=c("#1fbfc3", "#f57770"),
                       labels=c("Tropical", "Temperate")) +
  scale_x_continuous(breaks = 1:4+0.25, labels = c('Young', 'Moderate', 'Old', 'Unknown')) + 
  guides(color = guide_legend(override.aes = list(fill = "white"))) +
  theme(axis.title.y = element_markdown(size=10, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_um_group_count_plot2

# test age and UMR
LTR = LTR %>% mutate(UMR = case_when(umlen >= 0 ~ 'yes', TRUE ~ 'no'))
t.test(LTR$age ~ LTR$UMR)
t.test(subset(LTR, category == 'Tropical amplification')$age ~ subset(LTR, category == 'Tropical amplification')$UMR)

# UMR% fam stat
LTR_umr_sum2 = LTR %>%
  group_by(category) %>% summarise(category = unique(category),
                             um_rate = sum(umlen > 0, na.rm = T)/length(umlen),
                             age_mean = mean(age))

# LTR spec stat
table(LTR$category)
LTR %>% group_by(category) %>% 
  summarise(category = unique(category),
            ltrlen = median(ltrlen, na.rm = T),
            telen = median(ltrlen*2+intlen, na.rm = T),
            age_mean = median(age))

# stat LTR specs
LTR_ltrlen_cat_plot = ggplot(subset(LTR, category %in% c('Temperate removal', 'Tropical amplification')),
                             aes(category, ltrlen)) + geom_boxplot() + 
  scale_y_continuous(trans='log10') + theme_bw() + ylab('LTR region length (bp)') +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95))
LTR_ltrlen_cat_plot

LTR_telen_cat_plot = ggplot(subset(LTR, category %in% c('Temperate removal', 'Tropical amplification')),
                            aes(category, ltrlen*2+intlen)) + geom_boxplot() + 
  scale_y_continuous(trans='log10') + theme_bw() + ylab('LTR element length (bp)') +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95))
LTR_telen_cat_plot

LTR_domainct_cat_plot = ggplot(subset(subset(LTR, category %in% c('Temperate removal', 'Tropical amplification')),
                                      category %in% c('Temperate removal', 'Tropical amplification')),
                               aes(category, domain_ct)) + geom_boxplot() + 
  theme_bw() + ylab('Domain #') +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95))
LTR_domainct_cat_plot

table(LTR$domain_complete, LTR$category)
LTR_domain_cat_sum = data.frame(table(LTR$domain_complete, LTR$category)/
                                  rep(table(subset(LTR, domain_complete != 'NA')$category), each=2))
names(LTR_domain_cat_sum) = c('domain_complete', 'category', 'Frequency')
LTR_domain_cat_plot = ggplot(subset(LTR_domain_cat_sum, category %in% c('Temperate removal', 'Tropical amplification')),
                             aes(category, Frequency, fill = domain_complete)) + geom_col() +
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95))
LTR_domain_cat_plot

# add count to Fam_count
Fam_count$um_count = as.vector(table(subset(LTR, umlen>0)$id)[Fam_count$TE])
Fam_count$tot_count = as.vector(table(LTR_keep$id)[Fam_count$TE])
ggplot(subset(Fam_count, um_count>0), aes(tot_count, um_count)) + geom_point() #+ geom_line()
cor.test(Fam_count$um_count, log(Fam_count$tot_count))


####################################
#### genetic distance of LTRs  ####
####################################

# read composite recombination map
## Source data: https://datacommons.cyverse.org/browse/iplant/home/silastittes/parv_local_data/map/ogut_v5.map.txt
## This block of code was contributed by Silas Tittes to clean up and get cM/Mb from file
library(tidyverse)
recomb = read.table('ogut_v5.map.txt', sep = '\t', header = F)
gen_map_all_chr <- read_delim("ogut_v5.map.txt", delim = "\t") %>% 
  drop_na() %>%
  mutate(cm = cm + abs(min(cm))) %>%
  group_by(chr) %>% 
  group_modify(~{
    df1 <- slice(.x, -nrow(.x))
    df2 <- slice(.x, -1)
    to_keep <- df2$cm > df1$cm & df2$pos > df1$pos
    df1 <- df1[to_keep, ]
    df2 <- df2[to_keep, ]
    cm_mb <- tibble(cm_mb = 1e6*(df2$cm - df1$cm)/(df2$pos - df1$pos))
    bind_cols(df2, cm_mb)
  }) %>% 
  mutate(chr = paste0("chr", chr))
dim(recomb)
dim(gen_map_all_chr)

# get nearest cm/mb for each intact LTR
library(data.table)
new_LTR = data.frame()
for (i in 1:10){
  this_chr = subset(LTR, chr == paste('chr', i, sep = ''))
  this_recom = subset(gen_map_all_chr, chr == paste('chr', i, sep = ''))
  a=data.table(Value=this_chr$from)
  a[,merge:=Value]
  b=data.table(Value=this_recom$pos)
  b[,merge:=Value]
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  Merge_a_b = b[a,roll='nearest']
  Merge_a_b$cm_mb = this_recom$cm_mb[match(Merge_a_b$Value, this_recom$pos)]
  this_chr$cm_mb = Merge_a_b$cm_mb[match(this_chr$from, Merge_a_b$merge)]
  new_LTR = rbind(new_LTR, this_chr)
}
dim(new_LTR)
LTR$cm_mb = new_LTR$cm_mb[match(LTR$TE_id, new_LTR$TE_id)]
rm(new_LTR)

# get nearest cm/mb for each solo LTR
solo$chr = sub(".*_", '', solo$chr)
solo$condition = Fam_count$condition[match(solo$id, Fam_count$TE)]
solo$category = Fam_count$category[match(solo$id, Fam_count$TE)]
solo$status_class = Fam_count$status_class[match(solo$id, Fam_count$TE)]
solo$trop_temp_diff = Fam_count$trop_temp_diff[match(solo$id, Fam_count$TE)]
solo = solo %>% mutate(group = case_when(genome %in% NAM_trop ~ 'trop',
                                         genome %in% NAM_nontrop ~ 'temp',
                                         TRUE ~ 'admx'))
new_solo = data.frame()
for (i in 1:10){
  this_chr = subset(solo, chr == paste('chr', i, sep = ''))
  this_recom = subset(gen_map_all_chr, chr == paste('chr', i, sep = ''))
  a=data.table(Value=this_chr$from)
  a[,merge:=Value]
  b=data.table(Value=this_recom$pos)
  b[,merge:=Value]
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  Merge_a_b = b[a,roll='nearest']
  Merge_a_b$cm_mb = this_recom$cm_mb[match(Merge_a_b$Value, this_recom$pos)]
  this_chr$cm_mb = Merge_a_b$cm_mb[match(this_chr$from, Merge_a_b$merge)]
  new_solo = rbind(new_solo, this_chr)
}
dim(new_solo)
solo$cm_mb = new_solo$cm_mb[match(solo$id, new_solo$id)]
rm(new_solo)

# check data
p1 = ggplot(subset(LTR, genome == 'B97' & str_detect(chr, 'chr'))) + geom_point(aes(from, rrate)) + facet_wrap(~chr)
p2 = ggplot(gen_map_all_chr) + geom_point(aes(pos, cm_mb)) + facet_wrap(~chr) + ylim(0, 100)
ggplot(gen_map_all_chr) + geom_point(aes(pos, cm)) + facet_wrap(~chr)
p1 / p2

# stat group mean
tapply(LTR$cm_mb, LTR$condition, mean, na.rm = T)
tapply(subset(LTR, classification != 'LTR/unknown' & status_class != 'unknown')$cm_mb, 
       subset(LTR, classification != 'LTR/unknown' & status_class != 'unknown')$condition, mean, na.rm = T)

## permutation test on conditions
library('coin')
if(!require(rcompanion)){install.packages("rcompanion")}
# source for rcompanion: https://cran.microsoft.com/snapshot/2017-07-11/web/packages/rcompanion/index.html
install.packages("rcompanion")
if(!require(multcompView)){install.packages("multcompView")}

# pairwise test since it's significant
LTR$cm_mb.f = factor(LTR$cm_mb, ordered=TRUE)
LTR$trop_temp_diff = Fam_count$trop_temp_diff[match(LTR$id, Fam_count$TE)]

# construct a summary function
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm=TRUE)/sqrt(length(na.omit(x[[col]]))))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
  names(data_sum)[names(data_sum) == 'mean'] = varname
  return(data_sum)
}

# summarize mean and sd
LTR_cond_cm_mb_mean = groupwiseMean(cm_mb ~ category,
                                    data       = as.data.frame(subset(LTR, status_class != 'unknown' & 
                                                                        condition != 'other' & cm_mb > 0 &
                                                                        trop_temp_diff > 0)),
                                    conf       = 0.95,
                                    R          = 1000,
                                    percentile = TRUE,
                                    bca        = FALSE,
                                    digits     = 3)
LTR_cond_cm_mb_mean_plot = 
  ggplot(LTR_cond_cm_mb_mean, aes(x = category, y = Mean, color = category)) +
  geom_errorbar(aes(ymin = Percentile.lower, ymax = Percentile.upper), 
                width = 0.5, size  = 0.5) +
  geom_point(shape = 15, size  = 1, alpha = 2) + 
  scale_x_discrete(breaks = c('Temperate removal', 'Balanced (trop > temp)', 
                              'Tropical amplification', 'Drifting (trop > temp)'))  + 
  scale_color_manual(breaks = c('Temperate removal', 'Balanced (trop > temp)', 
                                'Tropical amplification', 'Drifting (trop > temp)'),
                     values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                     name = 'Classifications') + theme_bw() +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95),
        axis.title.y = element_text(size=8, face="bold"),
        plot.title = element_text(size=10),
        legend.position = "none") +
  labs(x='', y='Recombination rate (cM/Mb)', title = 'Intact LTR-RTs')
LTR_cond_cm_mb_mean_plot


# test solo recombination
solo_cond_cm_mb_mean = groupwiseMean(cm_mb ~ category,
                                    data       = as.data.frame(subset(solo, status_class != 'unknown' & 
                                                                        condition != 'other' & cm_mb > 0 &
                                                                        trop_temp_diff > 0)),
                                    conf       = 0.95,
                                    R          = 1000,
                                    percentile = TRUE,
                                    bca        = FALSE,
                                    digits     = 4)
solo_cond_cm_mb_mean_plot = 
  ggplot(solo_cond_cm_mb_mean, aes(x = category, y = Mean, color = category)) +
  geom_errorbar(aes(ymin = Percentile.lower, ymax = Percentile.upper), 
                width = 0.5, size  = 0.5) +
  geom_point(shape = 15, size  = 1, alpha = 2) + 
  scale_x_discrete(breaks = c('Temperate removal', 'Balanced (trop > temp)', 
                              'Tropical amplification', 'Drifting (trop > temp)'))  + 
  scale_color_manual(breaks = c('Temperate removal', 'Balanced (trop > temp)', 
                                'Tropical amplification', 'Drifting (trop > temp)'),
                     values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                     name = 'Classifications') + theme_bw() +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95),
        axis.title.y = element_text(size=8, face="bold"),
        plot.title = element_text(size=10),
        legend.position = "none") +
  labs(x='', y='Recombination rate (cM/Mb)', title = 'Solo LTRs')
solo_cond_cm_mb_mean_plot

# compare temp and trop
LTR_cond_cm_mb_group_mean = groupwiseMean(cm_mb ~ category + group,
                                           data       = as.data.frame(subset(LTR, cm_mb > 0 & group %in% c('temp', 'trop') &
                                                                               category %in% c('Temperate removal', 'Tropical amplification') &
                                                                               !(id %in% inconsistent_list))),
                                           conf       = 0.95,
                                           R          = 1000,
                                           percentile = TRUE,
                                           bca        = FALSE,
                                           digits     = 4)

LTR_cond_cm_mb_group_mean_plot = 
  ggplot(LTR_cond_cm_mb_group_mean, aes(x = category, y = Mean, color = group)) +
  geom_errorbar(aes(ymin = Percentile.lower, ymax = Percentile.upper), position = 'dodge',
                width = 0.5, size  = 0.5) + theme_bw() +
  geom_point(shape = 15, size  = 1, alpha = 2, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values=c('#53CACE', '#F5958F'), name = '', labels = c('Temp', 'Trop')) +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95),
        axis.title.y = element_text(size=8, face="bold"),
        plot.title = element_text(size=10),
        legend.position = "right") + 
  labs(x='', y='Recombination rate (cM/Mb)', title = 'Intact LTR-RTs')
LTR_cond_cm_mb_group_mean_plot

solo_cond_cm_mb_group_mean = groupwiseMean(cm_mb ~ category + group,
                                          data       = as.data.frame(subset(solo, cm_mb > 0 & group %in% c('temp', 'trop') &
                                                                              category %in% c('Temperate removal', 'Tropical amplification') &
                                                                              !(id %in% inconsistent_list))),
                                          conf       = 0.95,
                                          R          = 1000,
                                          percentile = TRUE,
                                          bca        = FALSE,
                                          digits     = 4)

solo_cond_cm_mb_group_mean_plot = 
  ggplot(solo_cond_cm_mb_group_mean, aes(x = category, y = Mean, color = group)) +
  geom_errorbar(aes(ymin = Percentile.lower, ymax = Percentile.upper), position = 'dodge',
                width = 0.5, size  = 0.5) + theme_bw() +
  geom_point(shape = 15, size  = 1, alpha = 2, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values=c('#53CACE', '#F5958F'), name = '', labels = c('Temp', 'Trop')) +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95),
        axis.title.y = element_text(size=8, face="bold"),
        plot.title = element_text(size=10),
        legend.position = "right") + 
  labs(x='', y='Recombination rate (cM/Mb)', title = 'Solo LTRs')
solo_cond_cm_mb_group_mean_plot


####################################
#### LTR distance to genes  ####
####################################
genedist = read.table('NAM.intact.LTR.genedist', sep = '\t')
names(genedist) = c('genome', 'chr', 'from', 'to', 'id', 'iden', 'geneL_chr', 
                    'geneL_from', 'geneL_to', 'geneL', 'geneL_dir', 'geneL_dist', 
                    'geneR_chr', 'geneR_from', 'geneR_to', 'geneR', 'geneR_dir', 'geneR_dist')
genedist$geneL_dist = as.numeric(genedist$geneL_dist)
genedist$geneR_dist = as.numeric(genedist$geneR_dist)
genedist$mindist = pmin(abs(genedist$geneL_dist), abs(genedist$geneR_dist))
genedist$TE_id = paste(genedist$chr, ":", genedist$from, "..", genedist$to, sep = '')
LTR$genedist = genedist$mindist[match(LTR$TE_id, genedist$TE_id)]
rm(genedist)

# different distances
options(scipen = 999)
LTR_genedist_hist = ggplot(LTR) + geom_histogram(aes(genedist+1)) + 
  scale_x_log10(breaks=c(1, 10, 1000, 100000), labels = c(0, 100, 10000, 1000000)) +
  labs(x='Distance to genes (bp)', y='All intact LTR-RTs (#)') + theme_bw()
LTR_genedist_hist
LTR_genedist_UMR_hist2 = ggplot(subset(LTR, umlen>0)) + geom_histogram(aes(genedist+1)) + 
  scale_x_log10(breaks=c(1, 10, 1000, 100000), labels = c(0, 100, 10000, 1000000)) +
  labs(x='Distance to genes (bp)', y='CHG UM intact LTR-RTs (#)') + theme_bw()
LTR_genedist_UMR_hist2


####################################
#### comparing trop temp diff.  ####
####################################
#stats
t.test(subset(LTR, group != 'admx')$iden ~ subset(LTR, group != 'admx')$group)
t.test(subset(LTR, group != 'admx')$age ~ subset(LTR, group != 'admx')$group)
t.test(subset(LTR, group != 'admx' & age < 20)$age ~ subset(LTR, group != 'admx' & age < 20)$group)
(0.9873941-0.9871096)/(2*3.3e-8)
t=11.301
2*pt(-abs(t),df=1126037) #exact p = 1.302102e-29
t.test(subset(LTR_nest_condense, group != 'admx')$nest_len ~ 
         subset(LTR_nest_condense, group != 'admx')$group) #p-value = 0.6458

# stats for fig. 4b
sum(subset(Fam_count, status_class == 'active')$trop_temp_diff) #26.06239 Mb
sum(subset(Fam_count, status_class == 'phasing down')$trop_temp_diff) #5.824402 Mb
sum(subset(Fam_count, status_class == 'inactive')$trop_temp_diff) #7.868869 Mb
sum(subset(Fam_count, status_class == 'active')$trop_temp_diff)/
  sum(subset(Fam_count, str_detect(Supfam, 'LTR'))$trop_temp_diff) #63.77809%
sum(subset(Fam_count, status_class == 'active')$trop_temp_diff)/
  sum(Fam_count$trop_temp_diff) #61.4727%
sum(subset(Fam_count, status_class == 'inactive')$trop_temp_diff)/
  sum(Fam_count$trop_temp_diff)#18.5601%
sum(subset(Fam_count, condition == 'Low remove High diff' & status_class == 'active' &
             Supfam == 'LTR/Gypsy')$trop_temp_diff) #14.94539 Mb
sum(subset(Fam_count, condition == 'Low remove High diff' & status_class == 'active' &
             Supfam == 'LTR/Gypsy')$trop_temp_diff)/
  sum(subset(Fam_count, str_detect(Supfam, 'LTR'))$trop_temp_diff) #36.57%
sum(subset(Fam_count, condition == 'High remove High diff' & status_class == 'active' &
             Supfam == 'LTR/Copia')$trop_temp_diff) #5.323094 Mb
sum(subset(Fam_count, condition == 'High remove High diff' & status_class == 'active' &
             Supfam == 'LTR/Copia')$trop_temp_diff)/
  sum(subset(Fam_count, str_detect(Supfam, 'LTR'))$trop_temp_diff) #13.03%
sum(subset(Fam_count, condition == 'Low remove High diff' & status_class == 'inactive' &
             Supfam == 'LTR/Gypsy')$trop_temp_diff) #4.023784 Mb

subset(Fam_count, condition == 'High remove High diff' & status_class == 'active' &
         +            Supfam == 'LTR/Copia' & trop_temp_diff > 2) #2.968885
subset(Fam_count, condition == 'Low remove High diff' & status_class == 'active' &
         Supfam == 'LTR/Gypsy' & trop_temp_diff > 2) #4.384978

histogram(subset(LTR, condition == 'Low remove High diff' & status_class == 'inactive' &
                   classification == 'LTR/Gypsy')$iden)
mean(subset(LTR, condition == 'Low remove High diff' & status_class == 'inactive' &
              classification == 'LTR/Gypsy')$iden) #0.9812941
(1-0.9812941)/(2*3.3e-8) #283422.7

# ridge plot of family age, divide trop and temp groups
LTR_mean_age = aggregate(sizelist_LTR$age, list(sizelist_LTR$id), mean)
LTR_mean_age$status_class = Fam_count$status_class[match(LTR_mean_age$Group.1, Fam_count$TE)]
LTR_mean_age$condition = Fam_count$condition[match(LTR_mean_age$Group.1, Fam_count$TE)]
LTR_mean_age = LTR_mean_age[with(LTR_mean_age, order(status_class, x)),]
sizelist_LTR$group = factor(sizelist_LTR$group, levels = c('trop', 'temp', 'admx'))
sizelist_LTR$id = factor(sizelist_LTR$id, levels = LTR_mean_age$Group.1) # sort ID with age
sizelist_LTR$category = LTR$category[match(sizelist_LTR$id, LTR$id)]

sizelist_LTR_ridge_plot_group2 = ggplot(sizelist_LTR[sizelist_LTR$group != 'admx',]) + 
  geom_density_ridges(aes(x = age, y = id, fill = group, height = ..density..),
                      stat = "density", scale = 10, alpha = 0.5, size = 0.25) +
  geom_brace(1450,1500,1,11, pointing="side", mid=0.5, lineend="round") + #curves for grouping
  geom_brace(1450,1500,12,22, pointing="side", mid=0.5, lineend="round") + #by ggbrace
  geom_brace(1450,1500,23,50, pointing="side", mid=0.5, lineend="round") +
  geom_point(data = LTR_mean_age, aes(x = 0, y = Group.1, color = condition), show.legend = F) +
  scale_fill_discrete(name="Group", breaks=c("trop", "temp"),
                      labels=c("Tropical", "Temperate")) +
  scale_color_manual(breaks=c("High remove High diff","High remove Low diff",
                              "Low remove High diff","Low remove Low diff"),
                     values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9')) +
  labs(x ="Insertion time (kya)", y = "Top 50 LTR families") + xlim(0, 1500) +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=6), 
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position=c(.8, .1),
        legend.key.size = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.3, 'cm'),
        plot.margin = margin(0.2, 1.5, 0.2, 0.2, "cm")
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
sizelist_LTR_ridge_plot_group2

# compare trop temp with consistent class
LTR_fam_info_group2$status_class = Fam_count$status_class[match(LTR_fam_info_group2$id, Fam_count$TE)]
LTR_fam_info_group_size_plot2 = ggplot(LTR_fam_info_group2, aes(status_class, size_mean, fill = group)) + 
  geom_boxplot(alpha = 0.5, size = 0.25) +
  scale_y_continuous(trans = 'log10') +
  scale_x_discrete(limits=c('active', 'phasing down', 'inactive'), labels = c("Young", 'Moderate', 'Old')) + 
  labs(x ="", y = "LTR family size (Mb)") +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95)) +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_fam_info_group_size_plot2

LTR_fam_info_group_count_plot2 = ggplot(LTR_fam_info_group2, aes(status_class, count, fill = group)) + 
  geom_boxplot(alpha = 0.5, size = 0.25) +
  scale_y_continuous(trans = 'log10') +
  scale_x_discrete(limits=c('active', 'phasing down', 'inactive'), labels = c("Young", 'Moderate', 'Old')) + 
  labs(x ="", y = "LTR family size (count)") +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95)) +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_fam_info_group_count_plot2

LTR_fam_info_group_age_plot2 = ggplot(LTR_fam_info_group2, aes(status_class, age_mean, fill = group)) + 
  geom_boxplot(alpha = 0.5, size = 0.25) +
  scale_y_continuous(trans = 'log10') +
  scale_x_discrete(limits=c('active', 'phasing down', 'inactive'), labels = c("Young", 'Moderate', 'Old')) + 
  labs(x ="", y = "LTR family age (kya)") +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95)) +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_fam_info_group_age_plot2


# summarize LTR fam of all genomes
LTR_count = subset(LTR, str_detect(id, ':', negate = T)) %>% group_by(id) %>%
  summarise(len_mean = mean(len, na.rm=T), age_mean = mean(age, na.rm=T),
            size_mean = sum(len, na.rm=T)/1000000/26)
Fam_count$age_mean = LTR_count$age_mean[match(Fam_count$TE, LTR_count$id)]

# draw SI - size, diff, and condition, differentiate trop larger or temp larger
LTR_fam_SI_size_tropbig_plot = ggplot(subset(Fam_count, category != 'Other' & trop_temp_diff > 0)) + 
  geom_point(aes(Size, SI, size = abs(trop_temp_diff), color = category)) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1, base = 2),
                     breaks=c(0, 2, 5, 10, 25, 70), limits = c(0, 80)) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 2),
                     breaks=c(0, 2, 10, 50, 200, 800), limits = c(0, 800)) +
  scale_color_manual(breaks = c("Temperate removal", "Balanced (trop > temp)",
                                "Tropical amplification", "Drifting (trop > temp)"),
                     values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                     name = 'Classifications') +
  scale_size(range = c(0, 10), limits = c(0, 5), name = 'Trop - Temp diff (Mb)') +
  labs(y ="Solo:intact ratio", x = "LTR family size (Mb)", title = 'Tropical larger families') +
  guides(color = guide_legend(order = 1), 
         size = guide_legend(order = 2, override.aes = list(color = "grey"))) +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_fam_SI_size_tropbig_plot

LTR_fam_SI_size_tempbig_plot = ggplot(subset(Fam_count, category != 'Other' & trop_temp_diff < 0)) + 
  geom_point(aes(Size, SI, size = abs(trop_temp_diff), color = category)) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1, base = 2),
                     breaks=c(0, 2, 5, 10, 25, 70), limits = c(0, 80)) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 2),
                     breaks=c(0, 2, 10, 50, 200, 800), limits = c(0, 800)) +
  scale_color_manual(breaks = c("Tropical removal", "Balanced (trop < temp)",
                                "Temperate amplification", "Drifting (trop < temp)"),
                     values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                     name = 'Classifications') +
  scale_size(range = c(0, 10), limits = c(0, 5), name = 'Temp - Trop diff (Mb)') +
  labs(y ="Solo:intact ratio", x = "LTR family size (Mb)", title = 'Temperate larger families') +
  guides(color = guide_legend(order = 1), 
         size = guide_legend(order = 2, override.aes = list(color = "grey"))) +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_fam_SI_size_tempbig_plot


# differentiate trop larger or temp larger, plot based on p-values
LTR_fam_SI_size_pval_tropbig_plot = ggplot(subset(Fam_count, category != 'Other' & trop_temp_diff > 0)) + 
  geom_point(aes(diff_pval, pval_si, size = abs(trop_temp_diff), color = category)) +
  geom_hline(yintercept=0.05, color = "black") +
  geom_vline(xintercept=0.05, color = "black") +
  scale_x_continuous(trans =  'log10', breaks=c(10^-6, 0.001, 0.05, 1), 
                     limits = c(10^-6, 1), labels=c(0, 0.001, 0.05, 1)) +
  scale_y_continuous(trans =  'log10', breaks=c(10^-6, 0.001, 0.05, 1), 
                     limits = c(10^-6, 1), labels=c(0, 0.001, 0.05, 1)) +
  scale_color_manual(breaks = c("Temperate removal", "Balanced (trop > temp)",
                                "Tropical amplification", "Drifting (trop > temp)"),
                     values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                     name = 'Classifications') +
  scale_size(range = c(0, 10), limits = c(0, 5), name = 'Trop - Temp diff (Mb)') +
  labs(y ="Solo:intact ratio (p value)", x = "LTR family size (p value)", title = 'Tropical larger families') +
  guides(color = guide_legend(order = 1), 
         size = guide_legend(order = 2, override.aes = list(color = "grey"))) +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_fam_SI_size_pval_tropbig_plot

LTR_fam_SI_size_pval_tempbig_plot = ggplot(subset(Fam_count, category != 'Other' & trop_temp_diff < 0)) + 
  geom_point(aes(diff_pval, pval_si, size = abs(trop_temp_diff), color = category)) +
  geom_hline(yintercept=0.05, color = "black") +
  geom_vline(xintercept=0.05, color = "black") +
  scale_x_continuous(trans =  'log10', breaks=c(10^-6, 0.001, 0.05, 1), 
                     limits = c(10^-6, 1), labels=c(0, 0.001, 0.05, 1)) +
  scale_y_continuous(trans =  'log10', breaks=c(10^-6, 0.001, 0.05, 1), 
                     limits = c(10^-6, 1), labels=c(0, 0.001, 0.05, 1)) +
  scale_color_manual(breaks = c("Tropical removal", "Balanced (trop < temp)",
                                "Temperate amplification", "Drifting (trop < temp)"),
                     values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                     name = 'Classifications') +
  scale_size(range = c(0, 10), limits = c(0, 5), name = 'Temp - Trop diff (Mb)') +
  labs(y ="Solo:intact ratio (p value)", x = "LTR family size (p value)", title = 'Temperate larger families') +
  guides(color = guide_legend(order = 1), 
         size = guide_legend(order = 2, override.aes = list(color = "grey"))) +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_fam_SI_size_pval_tempbig_plot


####################################
#### pan-genome LTR matrix  ####
####################################
## read data
synLTR = read.table('final_27_genome_TE_resolved_B73_added_v2.txt', header = T, sep = '\t')
#names(synLTR) = sub('_.*', '', na.omit(synLTR)[1,])

# correct genome ids
names(synLTR)[names(synLTR) == 'Oh7b'] <- 'Oh7B'
synLTR$Oh7B = sub('Oh7b', 'Oh7B', synLTR$Oh7B)

# remove Ab10
synLTR$B73AB10 = NULL

# remove admix
synLTR = select(synLTR, names(synLTR)[!names(synLTR) %in% NAM_adx])

# keep lines that contain intact or truncated LTRs 
count_NA_null <- function(x) {sum(str_count(x, 'null'), is.na(x), na.rm=T)}
synLTR = synLTR[apply(synLTR, 1, count_NA_null) != 23,]

# filter ;-containing rows (42 rows)
synLTR = synLTR[apply(synLTR, 1, function(x) {sum(str_count(x, ';'), na.rm=T)}) == 0, ]

# add info
synLTR$nonNA_count = apply(synLTR[,2:24], 1, function(x) {23-sum(is.na(x))})
synLTR$null_count = apply(synLTR[,2:24], 1, function(x) {sum(str_detect(x, 'null'), na.rm = T)})
synLTR$ltr_count = synLTR$nonNA_count - synLTR$null_count
synLTR$ins_freq = (synLTR$nonNA_count - synLTR$null_count)/synLTR$nonNA_count
synLTR$nonNA_trop = apply(synLTR[,NAM_trop], 1, function(x) {13-sum(is.na(x))})
synLTR$nonNA_temp = apply(synLTR[,NAM_nontrop], 1, function(x) {10-sum(is.na(x))})
synLTR$trop_count = synLTR$nonNA_trop - apply(synLTR[,NAM_trop], 1, function(x) {sum(str_detect(x, 'null'), na.rm = T)})
synLTR$temp_count = synLTR$nonNA_temp - apply(synLTR[,NAM_nontrop], 1, function(x) {sum(str_detect(x, 'null'), na.rm = T)})
synLTR$age = LTR$age[match(sub('.*_', '', sub( '_intact', '', apply(synLTR, 1, function(x) {str_subset(x, 'intact')[1]}))), LTR$TE_id)]
synLTR$genedist = LTR$genedist[match(sub('.*_', '', sub( '_intact', '', apply(synLTR, 1, function(x) {str_subset(x, 'intact')[1]}))), LTR$TE_id)]
synLTR$genetdist = LTR$genetdist[match(sub('.*_', '', sub( '_intact', '', apply(synLTR, 1, function(x) {str_subset(x, 'intact')[1]}))), LTR$TE_id)]
synLTR$rrate = LTR$rrate[match(sub('.*_', '', sub( '_intact', '', apply(synLTR, 1, function(x) {str_subset(x, 'intact')[1]}))), LTR$TE_id)]
synLTR$condition = LTR$condition[match(sub('.*_', '', sub( '_intact', '', apply(synLTR, 1, function(x) {str_subset(x, 'intact')[1]}))), LTR$TE_id)]
synLTR$category = LTR$category[match(sub('.*_', '', sub( '_intact', '', apply(synLTR, 1, function(x) {str_subset(x, 'intact')[1]}))), LTR$TE_id)]
synLTR$cm_mb = LTR$cm_mb[match(sub('.*_', '', sub( '_intact', '', apply(synLTR, 1, function(x) {str_subset(x, 'intact')[1]}))), LTR$TE_id)]
synLTR = synLTR %>%
  mutate(group = case_when(trop_count == 0 & temp_count > 0 ~ 'temp.uniq',
                           trop_count > 0 & temp_count == 0 ~ 'trop.uniq',
                           TRUE ~ 'shared'))
t.test(subset(synLTR, group == 'temp.uniq')$age, subset(synLTR, group == 'trop.uniq')$age) #188.6143 vs 185.0110, p-value < 2.2e-16
2*pt(-abs(8.5974),df=872539) # p = 8.167356e-18

## non-NA frequency
synLTR_nonNA_count_plot = ggplot(synLTR, aes(nonNA_count)) +
  geom_histogram(binwidth=1, center=0.5) +
  scale_y_log10() +
  labs(x ="Call frequency (n = 23)", y = "Number of syntenic LTR calls") +
  theme(legend.title=element_blank()) +
  theme_bw()
synLTR_nonNA_count_plot


## check NA effects
mincount = 2
sum(str_detect(unlist(subset(synLTR, nonNA_count >= mincount & nonNA_trop/13 >= mincount/23 & nonNA_temp/10 >= mincount/23)), 'intact'), na.rm=T)
synLTR_filterNA_SFS_plot1 = ggplot(subset(synLTR, nonNA_count >= mincount & nonNA_trop/13 >= mincount/23 & nonNA_temp/10 >= mincount/23), 
                                   aes(x = ins_freq, y = stat(count)/sum(count), fill = as.factor(nonNA_count))) + 
  geom_histogram(bins=0.1, binwidth = 0.05, position=position_dodge(0.05*0.7)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = 'LTR insertion frequency', y = 'Allele frequency', fill = "NonNA counts") + theme_bw()
synLTR_filterNA_SFS_plot1

mincount = 12
synLTR_filterNA_SFS_plot2 = ggplot(subset(synLTR, nonNA_count >= mincount & nonNA_trop/13 >= mincount/23 & nonNA_temp/10 >= mincount/23), 
                                   aes(x = ins_freq, y = stat(count)/sum(count), fill = as.factor(nonNA_count), color='black')) + 
  geom_histogram(bins=0.1, binwidth = 0.05) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values='black') +
  guides(color = FALSE, fill=guide_legend(ncol=2)) +
  labs(x = 'LTR insertion frequency', y = 'Allele frequency', fill = "NonNA counts") + theme_bw()
synLTR_filterNA_SFS_plot3 = ggplot(subset(synLTR, nonNA_count >= mincount & nonNA_trop/13 >= mincount/23 & nonNA_temp/10 >= mincount/23), 
                                   aes(x = ins_freq, y = stat(count)/sum(count), color = as.factor(nonNA_count))) + 
  geom_density() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  guides(color=guide_legend(ncol=2)) +
  labs(x = 'LTR insertion frequency', y = 'Population density', color = "NonNA counts") + theme_bw()
synLTR_filterNA_SFS_plot3


## filters using ≥50% call rate
synLTR_nonNA = subset(synLTR, nonNA_count == 23 & ins_freq > 0)
mincount = 12
synLTR_keep = subset(synLTR, nonNA_count >= mincount & nonNA_trop/13 >= mincount/23 & 
                       nonNA_temp/10 >= mincount/23 & ins_freq > 0)
synLTR_keep$fold_sfs = ifelse(synLTR_keep$ins_freq > 0.5, 1-synLTR_keep$ins_freq, synLTR_keep$ins_freq)


## freq and age for all synLTR elements
categories = seq(0, 1.0, by=0.1)
intact_freq <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(intact_freq) = c("TE_id", "ins_freq")
for (i in 1:(length(categories)-1)){
  this_list = str_subset(unlist(subset(synLTR_keep, ins_freq > categories[i] & 
                                         ins_freq <= categories[i+1])[,2:24]), 'intact')
  this_df = data.frame(TE_id = this_list, 
                       ins_freq = rep(paste('(', categories[i], ", ", categories[i+1], "]", sep=''), 
                                      length(this_list)))
  intact_freq = rbind(intact_freq, this_df)
}
intact_freq$age = LTR$age[match(sub('.*_', '', sub('_intact', '', intact_freq$TE_id)), LTR$TE_id)]
intact_freq_mean = tapply(intact_freq$age, intact_freq$ins_freq, mean, na.rm=T)
intact_freq_median = tapply(intact_freq$age, intact_freq$ins_freq, median, na.rm=T)
cor.test(as.numeric(factor(names(intact_freq_mean))), intact_freq_mean)
cor.test(as.numeric(factor(names(intact_freq_median))), intact_freq_median)

# stat age dist
synltr_freq_boxplot = ggplot(intact_freq[!is.na(intact_freq$age),], aes(ins_freq, age, group = ins_freq)) +
  geom_boxplot() +
  #geom_violin() +
  labs(x ="LTR insertion frequency", y = "Age of intact LTR-RTs (kya)") +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95),
        legend.title=element_blank())
synltr_freq_boxplot

## distance to genes
intact_freq$genedist = LTR$genedist[match(sub('.*_', '', sub( '_intact', '', intact_freq$TE_id)), LTR$TE_id)]
intact_freq$rrate = LTR$rrate[match(sub('.*_', '', sub( '_intact', '', intact_freq$TE_id)), LTR$TE_id)]
intact_freq$condition = LTR$condition[match(sub('.*_', '', sub( '_intact', '', intact_freq$TE_id)), LTR$TE_id)]
synLTR_genedist_boxplot = ggplot(na.omit(intact_freq), aes(ins_freq, (genedist+1)/1000, group = ins_freq)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw(base_size=12) +
  labs(x ="LTR insertion frequency", y = "Distance to nearest genes (kb)") +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95),
        legend.title=element_blank())
synLTR_genedist_boxplot
synLTR_genedist_mean = tapply(intact_freq$genedist, intact_freq$ins_freq, mean, na.rm = T)
synLTR_genedist_median = tapply(intact_freq$genedist, intact_freq$ins_freq, median, na.rm = T)
cor.test(synLTR_genedist_mean, as.numeric(order(names(synLTR_genedist_mean))))
cor.test(synLTR_genedist_median, as.numeric(order(names(synLTR_genedist_median))))
tapply(Fam_count$genedist_mean, Fam_count$condition, median, na.rm=T)
tapply(Fam_count$genedist_mean, Fam_count$status_class, median, na.rm=T)

## unique LTR for trop or temp
synLTR_keep = synLTR_keep %>%
  mutate(group = case_when(trop_count == 0 & temp_count > 0 ~ 'temp.uniq',
                           trop_count > 0 & temp_count == 0 ~ 'trop.uniq',
                           TRUE ~ 'shared'))
table(synLTR_keep$group)

# stat recombination
tapply(synLTR_keep$genetdist, synLTR_keep$group, mean, na.rm = T)
tapply(synLTR_keep$age, synLTR_keep$group, median, na.rm = T)
tapply(subset(synLTR_keep, group == 'temp.uniq')$genetdist, 
       subset(synLTR_keep, group == 'temp.uniq')$condition, mean, na.rm = T)
tapply(subset(synLTR_keep, group == 'trop.uniq')$genetdist, 
       subset(synLTR_keep, group == 'trop.uniq')$condition, mean, na.rm = T)

tapply(synLTR_keep$rrate, synLTR_keep$group, mean, na.rm = T)
tapply(subset(synLTR_keep, group == 'temp.uniq')$rrate, 
       subset(synLTR_keep, group == 'temp.uniq')$condition, mean, na.rm = T)
tapply(subset(synLTR_keep, group == 'trop.uniq')$rrate, 
       subset(synLTR_keep, group == 'trop.uniq')$condition, mean, na.rm = T)
tapply(subset(synLTR_keep, group == 'temp.uniq')$rrate, 
       subset(synLTR_keep, group == 'temp.uniq')$category, mean, na.rm = T)
tapply(subset(synLTR_keep, group == 'trop.uniq')$rrate, 
       subset(synLTR_keep, group == 'trop.uniq')$category, mean, na.rm = T)


t.test(synLTR_keep[synLTR_keep$group == 'temp.uniq', ]$age, synLTR_keep[synLTR_keep$group == 'trop.uniq', ]$age)
2*pt(-abs(9.2912),df=14593) #exact p val

# plot trop-temp uniq LTR age
t.test(subset(LTR, group != 'admx')$age ~ subset(LTR, group != 'admx')$group)
2*pt(-abs(-11.301),df=1126037) # p = 1.302102e-29

synLTR_uniq_boxplot_log = ggplot(synLTR_keep[!is.na(synLTR_keep$age),], aes(as.factor(group), age+1, fill = group)) +
  geom_boxplot(alpha = 0.5) + scale_y_continuous(trans = 'log10') +
  scale_fill_manual(name='', values = c("#999999", "#0072B2", "#D55E00")) +
  scale_x_discrete(labels = c('Shared', 'Temperate unique', 'Tropical unique')) +
  labs(x ="", y = "Age of intact LTR-RTs (kya)") +
  annotate("text", x = c(1, 2, 3), y = 2700,
           label = c('a', 'b', 'c')) +
  theme_bw(base_size=12) +
  theme(axis.title.y = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95),
        legend.position = "none")
synLTR_uniq_boxplot_log

mean(synLTR$nonNA_count)/23 #0.07677244
tapply(synLTR$age, synLTR$group, summary)
tapply(subset(synLTR, age<20)$age, subset(synLTR, age<20)$group, summary)
tapply(subset(synLTR, age<10)$age, subset(synLTR, age<10)$group, summary)
tapply(subset(synLTR, age<5)$age, subset(synLTR, age<5)$group, summary)
table(subset(synLTR, age<5)$group)
table(subset(synLTR, age<20)$group)
table(subset(synLTR, age>20)$group)
fisher.test(rbind(c(51944, 354189), c(74440, 461852)))
fisher.test(rbind(c(51944, 354189)/10, c(74440, 461852)/13)) #p-value = 4.532e-06, more uniq LTRs in trop that are <20kya

## trop and temp sfs
synLTR_keep$trop_freq = synLTR_keep$trop_count / synLTR_keep$nonNA_trop
synLTR_keep$temp_freq = synLTR_keep$temp_count / synLTR_keep$nonNA_temp

# test corr
cor.test(synLTR_keep$ins_freq, synLTR_keep$age)

## rare allele test
synLTR_rare_count = data.frame(trop = c(sum(subset(synLTR_keep, trop_freq > 0.2)$trop_count > 0),
                                        sum(subset(synLTR_keep, trop_freq <= 0.2)$trop_count > 0)),
                               temp = c(sum(subset(synLTR_keep, temp_freq > 0.2)$temp_count > 0),
                                        sum(subset(synLTR_keep, temp_freq <= 0.2)$temp_count > 0)))
rownames(synLTR_rare_count) = c('Allele freq ≥ 0.2', 'Allele freq < 0.2')


synLTR_rare_count = data.frame(trop = c(sum(subset(synLTR_keep, trop_freq > 0.8)$trop_count > 0),
                                        sum(subset(synLTR_keep, trop_freq > 0.2 & trop_freq <= 0.8)$trop_count > 0),
                                        sum(subset(synLTR_keep, trop_freq <= 0.2)$trop_count > 0)),
                               temp = c(sum(subset(synLTR_keep, temp_freq > 0.8)$temp_count > 0),
                                        sum(subset(synLTR_keep, temp_freq > 0.2 & temp_freq <= 0.8)$temp_count > 0),
                                        sum(subset(synLTR_keep, temp_freq <= 0.2)$temp_count > 0)))
rownames(synLTR_rare_count) = c('Allele freq (0.8, 1]', 'Allele freq (0.2, 0.8]', 'Allele freq (0, 0.2]')

synLTR_rare_count
fisher.test(synLTR_rare_count, workspace = 1200000)

# trop amp fams dist to genes
tapply(synLTR_trop_clean$genedist, synLTR_trop_clean$condition, mean, na.rm=T)
tapply(synLTR_temp_clean$genedist, synLTR_temp_clean$condition, mean, na.rm=T)


#################################################
#### Site Frequency Spectrum (SFS)  ####
#################################################
# output all
write.table(synLTR, file = "synLTR.txt", append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

# require nonNA rate >= 50%
write.table(synLTR_keep, file = "synLTR_keep.txt", append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

## test correlations
cor.test(synLTR_keep$ins_freq, synLTR_keep$age)
cor.test(synLTR_keep$ins_freq, synLTR_keep$genedist)
cor.test(synLTR_keep$ins_freq, synLTR_keep$genetdist)

## Non-NA allele frequency
synLTR_nonNA_SFS_plot = ggplot(synLTR_nonNA, aes(x = ins_freq, y = stat(count)/sum(count))) + 
  #  geom_histogram(bins=0.1, binwidth = 0.05) +
  geom_density() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(y = 'Population density', x = 'LTR frequency') + theme_bw()
synLTR_nonNA_SFS_plot

synLTR_keep_SFS_plot = ggplot(subset(synLTR_keep, ins_freq != 0), aes(x = ins_freq, y = stat(count)/sum(count))) + 
  #  geom_histogram(bins=0.1, binwidth = 0.05) +
  geom_density() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(y = 'Population density', x = 'LTR frequency') + theme_bw()
synLTR_keep_SFS_plot


## process SoFoS results
list_all = Sys.glob("./synLTR_keep.sfs/synLTR*sfs")
list_unfolded = list_all[grep("folded", list_all, invert = T)]
list_folded = list_all[grep("folded", list_all, invert = F)]
synLTR_keep.unfolded = lapply(list_unfolded, read.table, header = T, sep = ',', stringsAsFactors=F)
synLTR_keep.unfolded = data.frame(matrix(unlist(synLTR_keep.unfolded), ncol = length(synLTR_keep.unfolded) , byrow=F), stringsAsFactors=FALSE)
synLTR_keep.unfolded = synLTR_keep.unfolded[35:44,]
names(synLTR_keep.unfolded) = sub('./synLTR_keep.sfs/synLTR[._]', '' , list_unfolded) %>% sub('.sfs', '', .)
synLTR_keep.unfolded_sum = data.frame(rare = round(colSums(synLTR_keep.unfolded[1:2,])),
                                      common = round(colSums(synLTR_keep.unfolded[3:10,])))
synLTR_keep.unfolded = as.data.frame(sapply(synLTR_keep.unfolded, function(x){x/sum(x)}))
synLTR_keep.unfolded$number = as.numeric(rownames(synLTR_keep.unfolded))

synLTR_keep.folded = lapply(list_folded, read.table, header = T, sep = ',', stringsAsFactors=F)
synLTR_keep.folded = data.frame(matrix(unlist(synLTR_keep.folded), ncol = length(synLTR_keep.folded) , byrow=F), stringsAsFactors=FALSE)
synLTR_keep.folded = synLTR_keep.folded[20:24,]
names(synLTR_keep.folded) = sub('./synLTR_keep.sfs/synLTR[._]', '' , list_folded) %>% sub('.sfs', '', .)
synLTR_keep.folded_sum = data.frame(rare = round(colSums(synLTR_keep.folded[1,])),
                                    common = round(colSums(synLTR_keep.folded[2:5,])))
synLTR_keep.folded = as.data.frame(sapply(synLTR_keep.folded, function(x){x/sum(x)}))
synLTR_keep.folded$number = as.numeric(rownames(synLTR_keep.folded))


LTR_category_freq = data.frame(cbind(table(LTR$category)/sum(table(LTR$category)),
                                     table(synLTR$category)/sum(table(synLTR$category)),
                                     table(synLTR_keep$category)/sum(table(synLTR_keep$category)),
                                     table(subset(synLTR, !pan_TE %in% synLTR_keep$pan_TE)$category)/
                                       sum(table(subset(synLTR, !pan_TE %in% synLTR_keep$pan_TE)$category))))
names(LTR_category_freq) = c('All LTR', 'All synLTR', 'Filtered synLTR', 'Discarded synLTR')
LTR_category_freq$category = rownames(LTR_category_freq)



# plot out SFS group by temp and trop
synLTR_keep.sfs_group_mexnull_p = ggplot(synLTR_keep.unfolded %>% 
                                           select(number, keep_temp.mexnull, keep_trop.mexnull) %>%
                                           gather(., group, sfs, keep_temp.mexnull:keep_trop.mexnull, factor_key=TRUE)) + 
  geom_bar(aes(number, sfs, fill = group), stat = 'identity', position = 'identity', alpha = 0.4) + 
  labs(x="LTR insertion frequency", y="Proportion") +
  scale_x_continuous(breaks=seq(0, 10, by = 2)) + scale_fill_manual(values=c('#53CACE', '#F5958F'), labels = c('Temp', 'Trop'), name = '') +
  theme_bw()
synLTR_keep.sfs_group_mexnull_p

synLTR.convert.sfs_group_post20k_p = ggplot(synLTR_keep.unfolded %>% 
                                            select(number, convert.temp.post20k, convert.trop.post20k) %>%
                                            gather(., group, sfs, convert.temp.post20k:convert.trop.post20k, factor_key=TRUE)) + 
  geom_bar(aes(number, sfs, fill = group), stat = 'identity', position = 'identity', alpha = 0.4) + 
  labs(x="LTR insertion count", y="Proportion") +
  scale_x_continuous(breaks=seq(0, 10, by = 2)) + 
  scale_fill_manual(values=c('#53CACE', '#F5958F'), name = '', labels = c('Temp', 'Trop')) +
  theme_bw() +
  theme(axis.title.y = element_text(size=8, face="bold"),
        axis.title.x = element_text(size=8, face="bold"),
        legend.position = c(0.8, 0.8))
synLTR.convert.sfs_group_post20k_p

# for synLTR ≤50% missing
synLTR_keep.sfs_group_folded_p = ggplot(synLTR_keep.folded %>% 
                                          select(number, keep_temp.folded, keep_trop.folded) %>%
                                          gather(., group, sfs, keep_temp.folded:keep_trop.folded, factor_key=TRUE)) + 
  geom_bar(aes(number, sfs, fill = group), stat = 'identity', position = 'identity', alpha = 0.4) + 
  labs(x="Minor allele count", y="Proportion", title = 'LTR SFS') + ylim(c(0, 0.35)) +
  scale_x_continuous(breaks=seq(0, 5, by = 1)) + scale_fill_manual(values=c('#53CACE', '#F5958F'), labels=c('Temp', 'Trop'), name = '') +
  theme_bw()
synLTR_keep.sfs_group_folded_p

# profile SNP data
SNP_nonNA = read.table('merged_mexicana-v1_filtered-snps.biallelic.homo.troptemp.nonNA.count', header = F)
names(SNP_nonNA) = c('nonNA_count')

## process SNP SoFoS data
SNP_all = Sys.glob("./synLTR_keep.sfs/merged*sfs")
SNP_folded = SNP_all[grep("folded", SNP_all, invert = F)]

NAM.SNP.folded = lapply(SNP_folded, read.table, header = T, sep = ',', stringsAsFactors=F)
NAM.SNP.folded = data.frame(matrix(unlist(NAM.SNP.folded), ncol = length(NAM.SNP.folded) , byrow=F), stringsAsFactors=FALSE)
NAM.SNP.folded = NAM.SNP.folded[20:24,]
names(NAM.SNP.folded) = sub('./synLTR_keep.sfs/merged.*homo.', '' , SNP_folded) %>% sub('.folded.*', '', .)

NAM.SNP.folded_sum = data.frame(rare = round(colSums(NAM.SNP.folded[1,])),
                                    common = round(colSums(NAM.SNP.folded[2:5,])))
NAM.SNP.folded = as.data.frame(sapply(NAM.SNP.folded, function(x){x/sum(x)}))
NAM.SNP.folded$number = as.numeric(rownames(NAM.SNP.folded))

# plot
NAM.SNP.mexi_folded_p = ggplot(NAM.SNP.folded %>% 
                                  select(number, miss50.temp, miss50.trop) %>%
                                  gather(., group, sfs, miss50.temp:miss50.trop, factor_key=TRUE)) + 
  geom_bar(aes(number, sfs, fill = group), stat = 'identity', position = 'identity', alpha = 0.4) + 
  labs(x="Minor allele count", y="Proportion", title = 'SNP SFS (≤50% miss, mexi)') + ylim(c(0, 0.35)) +
  scale_x_continuous(breaks=seq(0, 5, by = 1)) + scale_fill_manual(values=c('#53CACE', '#F5958F'), labels=c('Temp', 'Trop'), name = '') +
  theme_bw()
NAM.SNP.mexi_folded_p


#################################################
#### family-based diff expression  ####
#################################################
load("~/Google Drive/study/Maize research/NAM/TE/expression/LTR_fam_tissue.RData")
str(LTR_fam_tissue)

write.table(LTR_fam_tissue, file = "LTR_fam_tissue_normalized.txt", append = FALSE, 
            quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = F, col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

# check overlaps
length(unique(LTR_fam_tissue$id))
length(unique(subset(LTR_fam_tissue, str_detect(Supfam, 'LTR/') & 
         status_class != 'unknown' & trop_temp_diff > 0)$id))
length(unique(subset(LTR, umlen>0 & status_class != 'unknown' & trop_temp_diff > 0)$id))

# all sig exp fams and UMR overlap
umr_ovlp_p = ggVennDiagram(list(`DE Fams` = unique(subset(LTR_fam_tissue, str_detect(Supfam, 'LTR/') & 
                                               status_class != 'unknown')$id), 
                   `UM Fams` = unique(subset(LTR, umlen>0 & 
                                         status_class != 'unknown')$id)),
              edge_size = 0.1) +
  scale_fill_gradient(low = "#E69F00", high = "#56B4E9") + 
  scale_color_grey(start = 0, end = 0) +
  theme(legend.position = 'none')
umr_ovlp_p

share_id = subset(unique(subset(LTR, umlen>0 & status_class != 'unknown')$id), 
                  unique(subset(LTR, umlen>0 & status_class != 'unknown')$id) %in% 
                    unique(subset(LTR_fam_tissue, str_detect(Supfam, 'LTR/') &  status_class != 'unknown')$id))
sum(subset(Fam_count, TE %in% share_id)$trop_temp_diff) #21.86349
sum(subset(Fam_count, str_detect(Supfam, 'LTR/'))$trop_temp_diff) #40.84318
21.86349/40.84318 #53.5%

# trop amp higher trop overlap with umr
umr_tropamp_ovlp_p = ggVennDiagram(list(`DE TropAmp. Fams` = unique(subset(LTR_fam_tissue, str_detect(Supfam, 'LTR/') &
                                                                             id %in% LTR_trop_big_only &
                                                                             trop_temp_diff > 0 & log2FoldChange > 0 &
                                                                             category == 'Tropical amplification')$id), 
                                        `UM TropAmp. Fams` = unique(subset(LTR, umlen>0 & genome %in% NAM_trop &
                                                                             category == 'Tropical amplification')$id)),
                                   edge_size = 0.1) +
  scale_fill_gradient(low = "#E69F00", high = "#56B4E9") + 
  scale_color_grey(start = 0, end = 0) +
  theme(legend.position = 'none')
umr_tropamp_ovlp_p

share_id = subset(unique(subset(LTR, umlen>0 & category == 'Tropical amplification')$id), 
                  unique(subset(LTR, umlen>0 & category == 'Tropical amplification')$id) %in% 
                    unique(subset(LTR_fam_tissue, str_detect(Supfam, 'LTR/') & log2FoldChange > 0 &
                                    category == 'Tropical amplification' &  trop_temp_diff > 0)$id))

trop_big = unique(subset(LTR_fam_tissue, padj<.05 & log2FoldChange>0)$id) #1948
temp_big = unique(subset(LTR_fam_tissue, padj<.05 & log2FoldChange<0)$id) #1620
trop_big_only = trop_big[!(trop_big %in% temp_big)] #1789
temp_big_only = temp_big[!(temp_big %in% trop_big)] #1461
LTR_trop_big_only = subset(trop_big_only, trop_big_only %in% Fam_count$TE[str_detect(Fam_count$Supfam, 'LTR/')])
LTR_temp_big_only = subset(temp_big_only, temp_big_only %in% Fam_count$TE[str_detect(Fam_count$Supfam, 'LTR/')])

# 15 tropical amplification families were consistnetly significantly more highly expressed in tropical genomes and possessed at least one UM- 5’-LTR
share_id = intersect(unique(subset(LTR, genome %in% NAM_trop & umlen>0 & category == 'Tropical amplification')$id),
                    unique(subset(LTR_fam_tissue, id %in% LTR_trop_big_only)$id))

sum(subset(Fam_count, TE %in% share_id)$trop_temp_diff) #11.47345
sum(subset(Fam_count, str_detect(Supfam, 'LTR/') & 
             category == 'Tropical amplification')$trop_temp_diff) #27.9494
13.46751/27.9494 #48.2%
13.46751/40.84318 #33.0%
12.62067/40.84318 #30.90031%
11.47345/40.84318 #28.09147%
sum(Fam_count$trop_temp_diff) #42.39669
sum(Fam_count$Size) #1796.746
sum(subset(Fam_count, (status_class == 'active' & 
             Supfam == 'LTR/Gypsy' & category == 'Tropical amplification') | 
             (status_class == 'active' & 
                Supfam == 'LTR/Copia' & category == 'Temperate removal'))$trop_temp_diff) #20.33408
20.33408/40.84318
20.33408/42.39669
sum(subset(Fam_count, (status_class == 'active' & 
                         Supfam == 'LTR/Gypsy' & category == 'Tropical amplification') | 
             (status_class == 'active' & 
                Supfam == 'LTR/Copia' & category == 'Temperate removal'))$Size) #232.1523
232.1523/1796.746 # = 12.92071%
sum(subset(Fam_count, status_class == 'active' & 
             Supfam %in% c('LTR/Gypsy','LTR/Copia'))$trop_temp_diff) #25.11275
25.11275/40.84318 #= 61.48579%
sum(subset(Fam_count, status_class == 'active' & 
             Supfam %in% c('LTR/Gypsy','LTR/Copia'))$Size) #25.11275
348.9117/1796.746 #= 19.41909%

LTR_fam_tissue$um_te_pcnt = Fam_count$um_te_pcnt[match(LTR_fam_tissue$id, Fam_count$TE)]
LTR_fam_tissue$um_diff = Fam_count$um_diff[match(LTR_fam_tissue$id, Fam_count$TE)]
cor.test(subset(LTR_fam_tissue, um_diff != 'NA' & um_diff != 'NaN' & um_diff != 'Inf')$um_diff,  #p-value = 0.04279, r= 0.1487106 
         subset(LTR_fam_tissue, um_diff != 'NA' & um_diff != 'NaN' & um_diff != 'Inf')$log2FoldChange)
cor.test(log2(subset(LTR_fam_tissue, um_diff != 0 & um_diff != 'NA' & um_diff != 'NaN' & um_diff != 'Inf')$um_diff),  #p-value = 0.06016, r= 0.151328 
         subset(LTR_fam_tissue, um_diff != 0 & um_diff != 'NA' & um_diff != 'NaN' & um_diff != 'Inf')$log2FoldChange)

LTR_tropbig_tissue_plot2 = ggplot(subset(Fam_count, trop_temp_diff > 0.025 &
                                           TE %in% LTR_trop_big_only), aes(Supfam, trop_temp_diff)) + 
  geom_bar(stat = 'identity', aes(fill = category, color = 'LTR family', size = 0.01), alpha = 1) + 
  scale_color_manual(values = "black", name = '') +
  scale_fill_manual(breaks = c('Temperate removal', 'Balanced (trop > temp)', 
                               'Tropical amplification', 'Drifting (trop > temp)'),
                    values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                    name = 'Classifications') +
  scale_size(range = c(0.1, 0.2), guide='none') +
  facet_grid(.~status_class, labeller = labeller(status_class = c('active' = "Young", 'phasing down' = 'Moderate', 'inactive' = 'Old'))) + 
  labs(x='') + ylab("Size diff for families sig. <br>higher expressed in Trop (Mb)") +
  guides(color = guide_legend(override.aes = list(fill = "white"))) +
  theme(axis.title.y = element_markdown(size=10, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_tropbig_tissue_plot2



####################################
#### gather everything together  ####
####################################
#### gather element-based info to big matrixes
str(LTR_sort)
str(LTR)
str(test_trop_temp)
str(LTR_fam_info)
LTR$status_class = LTR_fam_info$class[match(LTR$id, LTR_fam_info$id)]
LTR$domain_ct = LTR_sort$domain_ct[match(LTR$TE_id, LTR_sort$TE_id)]
LTR$domain_complete = LTR_sort$Complete[match(LTR$TE_id, LTR_sort$TE_id)]
LTR$sort_class = LTR_sort$Superfamily[match(LTR$TE_id, LTR_sort$TE_id)]

# convert ambiguous and sort
LTR_key <- list('Gypsy' = "LTR/Gypsy", 'Copia' = "LTR/Copia", 'unknown' = "LTR/unknown", 
                'mixture' = 'LTR/unknown', 'PiggyBac' = 'LTR/unknown', 'hAT' = 'LTR/unknown', 
                'MuDR_Mutator' = 'LTR/unknown')
LTR$sort_class = dplyr::recode(LTR$sort_class, !!!LTR_key)
LTR[is.na(LTR$sort_class),]$sort_class = 'LTR/unknown'
LTR[is.na(LTR$status_class),]$status_class = 'unknown'
LTR$status_class = factor(LTR$status_class, levels = c("active", "phasing down", "inactive", "unknown"))

#### gather family-based info to big matrixes
Fam_count$status_class = LTR_fam_info$class[match(Fam_count$TE, LTR_fam_info$id)]
Fam_count$temp_mean = test_trop_temp$temp_mean[match(Fam_count$TE, test_trop_temp$id)]
Fam_count$trop_mean = test_trop_temp$trop_mean[match(Fam_count$TE, test_trop_temp$id)]
Fam_count$diff_pval = test_trop_temp$pval[match(Fam_count$TE, test_trop_temp$id)]
Fam_count$trop_temp_diff = test_trop_temp$diff[match(Fam_count$TE, test_trop_temp$id)]
Fam_count$sort_class = LTR$sort_class[match(Fam_count$TE, LTR$id)]
Fam_count$sort_class <- as.factor(ifelse(is.na(Fam_count$sort_class), as.character(Fam_count$Supfam), Fam_count$sort_class))
Fam_count$sort_class = factor(Fam_count$sort_class, levels = c("LTR/Gypsy", "LTR/Copia", "LTR/unknown")) #order LTRs but make others NA

## make a new matrix
sum(trop_temp_diff$trop_mean_mb)
sum(subset(trop_temp_diff, str_detect(Supfam, 'LTR'))$trop_mean_mb)
sum(subset(test_trop_temp, str_detect(class, 'LTR'))$trop_mean)

trop_temp_diff_old = trop_temp_diff
trop_temp_diff = test_trop_temp[c('id', 'pval', 'temp_mean', 'trop_mean', 'diff')]
trop_temp_diff = merge(trop_temp_diff, 
                       trop_temp_si[c('id', 'pval', 'temp_mean', 'trop_mean', 'diff')],
                       by = c('id'), all.x = T)
trop_temp_diff = merge(trop_temp_diff,
                       trop_temp_fi[c('id', 'pval', 'temp_mean', 'trop_mean', 'diff')],
                       by = c('id'), all.x = T)
trop_temp_diff = merge(trop_temp_diff, 
                       trop_temp_age[c('id', 'pval', 'temp_mean', 'trop_mean', 'diff')],
                       by = c('id'), all.x = T)
colnames(trop_temp_diff) = c("id", "pval_mb", "temp_mean_mb", "trop_mean_mb", "diff_mb", 
                             "pval_si", "temp_mean_si", "trop_mean_si", "diff_si", 
                             "pval_fi", "temp_mean_fi", "trop_mean_fi", "diff_fi", 
                             "pval_kya", "temp_mean_kya", "trop_mean_kya", "diff_kya")
trop_temp_diff$status_class = LTR_fam_info$class[match(trop_temp_diff$id, LTR_fam_info$id)]
trop_temp_diff$Supfam = as.factor(class$Supfam[match(trop_temp_diff$id, sub('_LTR|_INT', '', class$TE))])
trop_temp_diff$fam_size = fam_size$mean[match(trop_temp_diff$id, fam_size$TE)]
str(trop_temp_diff)

## convert NAs and order factors
Fam_count[is.na(Fam_count$status_class),]$status_class = 'unknown'
Fam_count$status_class = factor(Fam_count$status_class, levels = c("active", "phasing down", "inactive", "unknown"))
trop_temp_diff[is.na(trop_temp_diff$status_class),]$status_class = 'unknown'
trop_temp_diff$Supfam = factor(trop_temp_diff$Supfam, levels = c("LTR/Gypsy", "LTR/Copia", "LTR/unknown"))
trop_temp_diff[is.na(trop_temp_diff$Supfam),]$Supfam = 'LTR/unknown'
trop_temp_diff$status_class = factor(trop_temp_diff$status_class, levels = c("active", "phasing down", "inactive", "unknown"))


## add test judgements
trop_temp_diff = trop_temp_diff %>%
  mutate(condition = case_when(pval_si <= 0.05 & pval_mb <= 0.05 ~ 'High remove High diff',
                               pval_si > 0.05 & pval_mb <= 0.05 ~ 'Low remove High diff',
                               pval_si <= 0.05 & pval_mb > 0.05 ~ 'High remove Low diff',
                               pval_si > 0.05 & pval_mb > 0.05 ~ 'Low remove Low diff',
                               TRUE ~ 'Other'))
trop_temp_diff = trop_temp_diff %>%
  mutate(category = case_when(condition == 'High remove High diff' & diff_mb > 0 ~ 'Temperate removal',
                              condition == 'Low remove High diff' & diff_mb > 0 ~  'Tropical amplification',
                              condition == 'High remove Low diff' & diff_mb > 0 ~ 'Balanced (trop > temp)',
                              condition == 'Low remove Low diff' & diff_mb > 0 ~ 'Drifting (trop > temp)',
                              condition == 'High remove High diff' & diff_mb < 0 ~ 'Tropical removal',
                              condition == 'Low remove High diff' & diff_mb < 0 ~  'Temperate amplification',
                              condition == 'High remove Low diff' & diff_mb < 0 ~ 'Balanced (trop < temp)',
                              condition == 'Low remove Low diff' & diff_mb < 0 ~ 'Drifting (trop < temp)',
                              TRUE ~ 'Other'))
Fam_count$condition = trop_temp_diff$condition[match(Fam_count$TE, trop_temp_diff$id)]
Fam_count[is.na(Fam_count$condition),]$condition = 'other'
LTR$condition = Fam_count$condition[match(LTR$id, Fam_count$TE)]
LTR[is.na(LTR$condition),]$condition = 'other'
Fam_count$category = trop_temp_diff$category[match(Fam_count$TE, trop_temp_diff$id)]
Fam_count[is.na(Fam_count$category),]$category = 'other'
LTR$category = Fam_count$category[match(LTR$id, Fam_count$TE)]
LTR[is.na(LTR$category),]$category = 'other'
Fam_count$pval_si = trop_temp_diff$pval_si[match(Fam_count$TE, trop_temp_diff$id)]
sum(subset(trop_temp_diff, condition=='High remove High diff' & diff_si<0)$diff_mb) #7.3234 Mb
sum(subset(trop_temp_diff, category=='Temperate removal')$diff_mb) #8.000 Mb

# trop rmv < trop amp, trop size > temp
subset(trop_temp_diff, pval_si<0.05 & diff_si >0 & pval_mb < 0.05 & diff_mb > 0)
sum(subset(trop_temp_diff, pval_si<0.05 & diff_si >0 & pval_mb < 0.05 & diff_mb > 0)$diff_mb)

# temp rmv < temp amp, temp size > temp
table(trop_temp_diff$category)
sum(subset(trop_temp_diff, category == 'Tropical removal' & diff_si < 0)$diff_mb) #incon, bad, n = 11
sum(subset(trop_temp_diff, category == 'Tropical removal' & diff_si > 0)$diff_mb)
sum(subset(trop_temp_diff, category == 'Temperate removal' & diff_si < 0)$diff_mb) 
sum(subset(trop_temp_diff, category == 'Temperate removal' & diff_si > 0)$diff_mb) #incon, bad, n = 20
sum(subset(trop_temp_diff, category == 'Temperate amplification' & diff_si < 0)$diff_mb) #incon, ok
sum(subset(trop_temp_diff, category == 'Temperate amplification' & diff_si > 0)$diff_mb)
sum(subset(trop_temp_diff, category == 'Tropical amplification' & diff_si < 0)$diff_mb) 
sum(subset(trop_temp_diff, category == 'Tropical amplification' & diff_si > 0)$diff_mb) #incon, ok
sum(subset(trop_temp_diff, category == 'Balanced (trop < temp)' & diff_si < 0)$diff_mb) #incon, bad
sum(subset(trop_temp_diff, category == 'Balanced (trop < temp)' & diff_si > 0)$diff_mb)
sum(subset(trop_temp_diff, category == 'Balanced (trop > temp)' & diff_si < 0)$diff_mb) 
sum(subset(trop_temp_diff, category == 'Balanced (trop > temp)' & diff_si > 0)$diff_mb) #incon, bad
sum(subset(trop_temp_diff, category == 'Drifting (trop < temp)' & diff_si < 0)$diff_mb) #incon, ok
sum(subset(trop_temp_diff, category == 'Drifting (trop < temp)' & diff_si > 0)$diff_mb)
sum(subset(trop_temp_diff, category == 'Drifting (trop > temp)' & diff_si < 0)$diff_mb) 
sum(subset(trop_temp_diff, category == 'Drifting (trop > temp)' & diff_si > 0)$diff_mb) #incon, ok

# all inconsistent categories
subset(trop_temp_diff, (diff_si > 0 & diff_mb < 0) | (diff_si < 0 & diff_mb > 0))
dim(subset(trop_temp_diff, diff_si*diff_mb < 0))
sum(subset(trop_temp_diff, diff_si*diff_mb < 0)$diff_mb)
sum(subset(trop_temp_diff, diff_si*diff_mb > 0)$diff_mb)
dim(trop_temp_diff)
table((subset(trop_temp_diff, diff_si*diff_mb < 0)$category))
sum(subset(trop_temp_diff, diff_si*diff_mb < 0 & category == 'Temperate removal')$diff_mb)
subset(trop_temp_diff, diff_si*diff_mb < 0 & category == 'Temperate removal')

# stats
ggplot(Fam_count, aes(status_class, trop_temp_diff)) + geom_point()
tapply(subset(LTR_fam_info_group, class == 'active')$size_mean, 
       subset(LTR_fam_info_group, class == 'active')$group, sum)
tapply(subset(LTR_fam_info_group, class == 'phasing down')$size_mean, 
       subset(LTR_fam_info_group, class == 'phasing down')$group, sum)
tapply(subset(LTR_fam_info_group, class == 'inactive')$size_mean, 
       subset(LTR_fam_info_group, class == 'inactive')$group, sum)
tapply(LTR_fam_info_group$size_mean, LTR_fam_info_group$group, sum)
tapply(trop_temp_diff$fam_size, trop_temp_diff$status_class, sum)
tapply(trop_temp_diff$temp_mean_mb, trop_temp_diff$status_class, sum)
tapply(trop_temp_diff$trop_mean_mb, trop_temp_diff$status_class, sum)
tapply(trop_temp_diff$diff_mb, trop_temp_diff$status_class, sum)
tapply(trop_temp_diff$diff_mb, trop_temp_diff$condition, sum)
tapply(trop_temp_diff$diff_mb, trop_temp_diff$category, sum)
table(trop_temp_diff$category)
sum(subset(trop_temp_diff, !id %in% LTR_fam_info_group$id)$diff_mb)
sum(subset(trop_temp_diff, !id %in% LTR_fam_info_group$id & diff_mb < 0)$diff_mb)
sum(subset(trop_temp_diff, !id %in% LTR_fam_info_group$id & diff_mb > 0)$diff_mb)
tapply(subset(trop_temp_diff, Supfam == 'LTR/Gypsy')$temp_mean_si, subset(trop_temp_diff, Supfam == 'LTR/Gypsy')$status_class, sum)/tapply(subset(trop_temp_diff, Supfam == 'LTR/Copia')$temp_mean_si, subset(trop_temp_diff, Supfam == 'LTR/Copia')$status_class, sum)
tapply(subset(trop_temp_diff, Supfam == 'LTR/Gypsy')$temp_mean_si, subset(trop_temp_diff, Supfam == 'LTR/Gypsy')$status_class, mean)/tapply(subset(trop_temp_diff, Supfam == 'LTR/Copia')$temp_mean_si, subset(trop_temp_diff, Supfam == 'LTR/Copia')$status_class, mean)
tapply(trop_temp_diff$diff_mb, trop_temp_diff$condition, sum)
dim(subset(Fam_count, Size < 0.1))
dim(subset(Fam_count, Size < 0.1 & str_detect(Supfam, 'LTR')))
dim(subset(LTR, id %in% subset(Fam_count, Size < 0.1)$TE))
length(unique(subset(LTR, id %in% subset(Fam_count, Size < 0.1)$TE)$id))
sum(subset(trop_temp_diff, Supfam == 'LTR/Gypsy' & status_class == 'active' & category == 'Tropical amplification')$diff_mb)
sum(subset(trop_temp_diff, Supfam == 'LTR/Copia' & status_class == 'active' & category == 'Temperate removal')$diff_mb)
sum(subset(trop_temp_diff, str_detect(Supfam, 'LTR'))$diff_mb)
sum(subset(trop_temp_diff, status_class == 'active')$diff_mb)
14.99266/42.36213 #35.4%
5.341415/42.36213 #12.6%
26.06267/42.36213 #61.5%
sum(Fam_count$trop_temp_diff_int, na.rm = T) #19.60124
sum(subset(Fam_count, status_class == 'active')$trop_temp_diff_int, na.rm = T) #18.7726
18.7726/19.60124
test = subset(trop_temp_diff, category == 'Tropical removal' & diff_si > 0)$id # 2 fams
test = subset(trop_temp_diff, category == 'Temperate removal' & diff_si < 0)$id # 16
test = subset(trop_temp_diff, category == 'Temperate amplification' & diff_si > 0)$id # 25 fams
test = subset(trop_temp_diff, category == 'Tropical amplification' & diff_si < 0)$id # 62 fams
dim(subset(LTR, id %in% test))
table(subset(LTR, id %in% test)$id)
table(subset(LTR, id == 'vegu_AC190718_85')$domain_ct)
mean(subset(LTR, id == 'vegu_AC190718_85')$age) #284kya (overall mean = 192kya)

### make plots
# total size
LTR_status_diff_plot2 = ggplot(subset(Fam_count, str_detect(Supfam, 'LTR/') & status_class != 'unknown'), aes(Supfam, trop_temp_diff)) + 
  #LTR_status_diff_plot2 = ggplot(Fam_count[str_detect(Fam_count$Supfam, 'LTR/'),], aes(sort_class, trop_temp_diff)) + 
  geom_bar(stat = 'identity', aes(fill = condition, color = 'LTR family', size = 0.01), alpha = 1) + 
  scale_color_manual(values = "black", name = '') +
  guides(color = guide_legend(override.aes = list(fill = "white"))) +
  scale_fill_manual(labels = c('Removal', 'Trop-Temp balanced', 
                               'Amplification', 'Random drifting'),
                    values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                    name = 'Classifications') +
  scale_size(range = c(0.001, 0.15), guide='none') +
  facet_grid(.~status_class, labeller = labeller(status_class = c('active' = "Young", 'phasing down' = 'Moderate', 'inactive' = 'Old'))) + 
  labs(x='', y = "Trop - Temp Fam Size (Mb)") +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_status_diff_plot2

intLTR_status_diff_plot2 = ggplot(subset(Fam_count, str_detect(Supfam, 'LTR/') & status_class != 'unknown'), aes(Supfam, trop_temp_diff_int)) + 
  geom_bar(stat = 'identity', aes(fill = condition, color = 'LTR family', size = 0.01), alpha = 1) + 
  scale_color_manual(values = "black", name = '') +
  guides(color = guide_legend(override.aes = list(fill = "white"))) +
  scale_fill_manual(labels = c('Removal', 'Trop-Temp balanced', 
                               'Amplification', 'Random drifting'),
                    values = c('#FABEBB', '#D2D296', '#9DDEC2', '#9BD9F9'),
                    name = 'Classifications') +
  scale_size(range = c(0.001, 0.15), guide='none') +
  facet_grid(.~status_class, labeller = labeller(status_class = c('active' = "Young", 'phasing down' = 'Moderate', 'inactive' = 'Old'))) + 
  labs(x='', y = "Trop - Temp Intact Size (Mb)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=12, angle=35, vjust=1, hjust=0.95)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
intLTR_status_diff_plot2

## test if differential size families are enrihed in the young group
sum(Fam_count$trop_temp_diff_int, na.rm=T) #19.58237
sum(subset(Fam_count, int_diff_pval<0.05)$trop_temp_diff_int, na.rm=T) #18.10903
sum(subset(Fam_count, status_class == 'active')$trop_temp_diff_int, na.rm=T) #18.76262
18.76262/19.58237 = 95.81384%
sum(subset(LTR, status_class == 'active')$len)/sum(LTR$len) #4576427049/13645407758 = 33.53822%

diff_size_stat = as.data.frame(cbind(table(LTR$status_class),
                                     table(subset(LTR, id %in% subset(Fam_count, int_diff_pval<0.05)$TE)$status_class)))
names(diff_size_stat) = c('All', 'DS_fams')
diff_size_stat$status_class = rownames(diff_size_stat)
fisher.test(rbind(c(436332, 87370), c(490762, 809103))) #p-value < 2.2e-16

diff_size_stat_plot = ggplot(gather(subset(diff_size_stat, status_class != 'unknown'), Group, Count, All:DS_fams, factor_key=TRUE), 
       aes(status_class, Count, fill=Group)) +
  scale_x_discrete(labels = c("Young", 'Moderate', 'Old')) +
  geom_col(position = 'dodge') + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=12))
diff_size_stat_plot


# plot stacked values
# make a new df trop_temp_diff_size
trop_temp_diff_size = reshape2::melt(trop_temp_diff, id.vars='id', variable.name = 'group', value.name = 'fam_size') #lots of columns merged
trop_temp_diff_size = trop_temp_diff_size[str_detect(trop_temp_diff_size$group, 'mean_mb'),] #filter for the rows containing fam size
trop_temp_diff_size$group = sub('_mean_mb', '', trop_temp_diff_size$group)
trop_temp_diff_size$Supfam = trop_temp_diff$Supfam[match(trop_temp_diff_size$id, trop_temp_diff$id)]
trop_temp_diff_size$fam_size = as.numeric(trop_temp_diff_size$fam_size)
trop_temp_diff_size$status_class = LTR_fam_info_group$class[match(paste(trop_temp_diff_size$id, trop_temp_diff_size$group, sep = '_'),
                                                                  LTR_fam_info_group$id_group)]
trop_temp_diff_size$id_group = paste(trop_temp_diff_size$id, trop_temp_diff_size$group, sep = "_")
LTR_fam_info_group$id_group = paste(LTR_fam_info_group$id, LTR_fam_info_group$group, sep = "_")
trop_temp_diff_size$status_class = LTR_fam_info_group$class[match(trop_temp_diff_size$id_group, LTR_fam_info_group$id_group)]
table(LTR_fam_info_group$class, LTR_fam_info_group$group, LTR_fam_info_group$Supfam)
table(trop_temp_diff_size$status_class, trop_temp_diff_size$group, trop_temp_diff_size$Supfam)
tapply(subset(trop_temp_diff, Supfam == 'LTR/Copia')$temp_mean_si, subset(trop_temp_diff, Supfam == 'LTR/Copia')$status_class, sum)
tapply(subset(trop_temp_diff, Supfam == 'LTR/Gypsy')$temp_mean_si, subset(trop_temp_diff, Supfam == 'LTR/Gypsy')$status_class, sum)

LTR_status_group_si_plot = ggplot() + 
  geom_bar(data = trop_temp_diff[trop_temp_diff$id %in% LTR_count_all[LTR_count_all$id_count>=10,]$id,], 
           stat = 'identity', alpha = 0.5, width = 0.35, 
           aes(as.numeric(Supfam)+0.05, trop_mean_si, fill = '#1fbfc3', color = 'LTR family', size = 0.05)) +
  geom_bar(data = trop_temp_diff[trop_temp_diff$id %in% LTR_count_all[LTR_count_all$id_count>=10,]$id,],
           stat = 'identity', alpha = 0.5, width = 0.35, 
           aes(as.numeric(Supfam)+0.45, temp_mean_si, fill = '#f57770', color = 'LTR family', size = 0.05)) +
  labs(x ="", y = "Cumulative solo:intact ratio") +
  scale_size(range = c(0.0001, 0.15), guide='none') +
  guides(color = guide_legend(override.aes = list(fill = "white"))) +
  facet_grid(.~status_class, labeller = labeller(status_class = c('active' = "Young", 'phasing down' = 'Moderate', 'inactive' = 'Old'))) + 
  scale_color_manual(values = "black", name = '') +
  scale_x_continuous(breaks = 1:3+0.25, labels = unique(na.omit(LTR_fam_info_group$Supfam))) + 
  scale_fill_discrete(name="Group", breaks=c("#1fbfc3", "#f57770"),
                      labels=c("Tropical", "Temperate"),
                      guide = guide_legend(reverse=TRUE)) +
  theme(axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
LTR_status_group_si_plot

# summary
solo_intact_mean = data.frame(id = solo_intact$id, mean_si = rowMeans(solo_intact[2:26]))
Fam_count$mean_si = solo_intact_mean$mean_si[match(Fam_count$TE, solo_intact_mean$id)]

tapply(trop_temp_diff$diff_mb, trop_temp_diff$condition, function(x)sum(x))
sum(Fam_count[Fam_count$condition == 'other',]$trop_temp_diff)
length(Fam_count[Fam_count$condition == 'other',]$trop_temp_diff)
tapply(trop_temp_diff[trop_temp_diff$condition == 'High remove High diff',]$diff_mb, trop_temp_diff[trop_temp_diff$condition == 'High remove High diff',]$Supfam, function(x)sum(x))
tapply(trop_temp_diff[trop_temp_diff$condition == 'High remove High diff',]$diff_mb, trop_temp_diff[trop_temp_diff$condition == 'High remove High diff',]$status_class, function(x)sum(x))
tapply(trop_temp_diff[trop_temp_diff$condition == 'Low remove High diff',]$diff_mb, trop_temp_diff[trop_temp_diff$condition == 'Low remove High diff',]$Supfam, function(x)sum(x))
tapply(trop_temp_diff[trop_temp_diff$condition == 'Low remove High diff',]$diff_mb, trop_temp_diff[trop_temp_diff$condition == 'Low remove High diff',]$status_class, function(x)sum(x))
t.test(LTR[LTR$id %in% trop_temp_diff[trop_temp_diff$condition == 'High remove High diff',]$id,]$iden,
       LTR[LTR$id %in% trop_temp_diff[trop_temp_diff$condition == 'Low remove High diff',]$id,]$iden)
tapply(subset(trop_temp_diff_size, status_class == 'active')$fam_size, 
       subset(trop_temp_diff_size, status_class == 'active')$group, sum)
tapply(subset(trop_temp_diff_size, status_class == 'phasing down')$fam_size, 
       subset(trop_temp_diff_size, status_class == 'phasing down')$group, sum)
tapply(subset(trop_temp_diff_size, status_class == 'inactive')$fam_size, 
       subset(trop_temp_diff_size, status_class == 'inactive')$group, sum)
tapply(trop_temp_diff_size$fam_size, trop_temp_diff_size$group, sum, na.rm = T)
tapply(subset(trop_temp_diff_size, status_class == 'phasing down' & Supfam == 'LTR/Gypsy')$fam_size, 
       subset(trop_temp_diff_size, status_class == 'phasing down' & Supfam == 'LTR/Gypsy')$group, sum)
tapply(subset(trop_temp_diff_size, status_class == 'phasing down' & Supfam == 'LTR/Copia')$fam_size, 
       subset(trop_temp_diff_size, status_class == 'phasing down' & Supfam == 'LTR/Copia')$group, sum)


####################################
#### Paper plots  ####
####################################
## Suppl. Fig. 1
suppFig1 = rare_LTR_plot2 + Fam_size_freq_int_plot_log + Fam_size_freq_homo_plot_log +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 12, face = 'bold'))
# plot out
pdf("suppFig 1.pdf", width=10,height=9,pointsize=12, paper='special')
suppFig1
dev.off()

pdf("suppFig 2.pdf", width=10,height=4,pointsize=12, paper='special')
fam_size_rank_plot
dev.off()

pdf("suppFig 3.pdf", width=10,height=4,pointsize=12, paper='special')
Fam_count_freq_plot
dev.off()

pdf("suppFig 4.pdf", width=15,height=6,pointsize=12, paper='special')
snp_PCA_biplot_group + busco_tree_p + 
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 18, face = 'bold'))
dev.off()

pdf("suppFig 5.pdf", width=12,height=7,pointsize=12, paper='special')
trop_temp_size_int_cmb_plot + intLTR_status_diff_plot2 + 
  diff_size_stat_plot + plot_spacer() +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 18, face = 'bold'))
dev.off()

## suppFig 6
suppFig6 = LTR_fam_SI_size_tropbig_plot + LTR_fam_SI_size_tempbig_plot + 
  LTR_status_group_si_plot + LTR_fam_info_group_age_plot2 + 
  plot_layout(design = 'AB\nCD', heights = c(3, 2)) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold'))
# plot out
pdf("suppFig6.pdf", width=12,height=7,pointsize=12, paper='special')
suppFig6
dev.off()

## suppFig 7
suppFig7 = (LTR_ltrlen_cat_plot | LTR_telen_cat_plot | 
               LTR_domainct_cat_plot | LTR_domain_cat_plot) + plot_layout(ncol = 4) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))
# plot out
pdf("suppFig 7.pdf", width=10,height=3,pointsize=12, paper='special')
suppFig7
dev.off()

# suppl fig 8
supplFig8 = (UM_LTR_topD_p1 + ggtitle('CRM2_7577nt') | UM_LTR_topD_p2 | UM_LTR_topD_p3) / 
  (UM_LTR_topC_p1 + ggtitle('uwum_AC213069_12092') | UM_LTR_topC_p2 | UM_LTR_topC_p3)  /
  (UM_LTR_topA_p1 + ggtitle('xilon_diguus_AC203313_7774') | UM_LTR_topA_p2 | UM_LTR_topA_p3) /
  (UM_LTR_topB_p1 + ggtitle('ji_AC215728_13156') | UM_LTR_topB_p2 | UM_LTR_topB_p3) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))
# plot out
pdf("supplFig8.pdf", width=10,height=11,pointsize=12, paper='special')
supplFig8
dev.off()

## suppFig 9
layout <- "
AB
CC"
suppFig9 = LTR_genedist_hist + LTR_genedist_UMR_hist2 +
  LTR_umlen_plot +
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold'))
# plot out
pdf("suppFig9.pdf", width=8,height=6,pointsize=12, paper='special')
suppFig9
dev.off()

## suppFig 12
suppFig12 = (umr_ovlp_p | umr_tropamp_ovlp_p) + 
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))
# plot out
pdf("suppFig12_v2.pdf", width=12,height=3,pointsize=12, paper='special')
suppFig12
dev.off()

## suppFig 14
layout <- "
AA
BB
CC
DE"
suppFig14 = (synLTR_nonNA_count_plot + synLTR_keep_SFS_plot + synLTR_keep_SFS_post20k_plot) - 
  synLTR_filterNA_SFS_plot1 +
  (synLTR_filterNA_SFS_plot2 + theme(legend.position = "none") + synLTR_filterNA_SFS_plot3) +
  synltr_freq_boxplot + synLTR_genedist_boxplot + plot_layout(design = layout) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))
# plot out
pdf("suppFig14.pdf", width=10,height=11,pointsize=12, paper='special')
suppFig14
dev.off()


# suppl fig 15
supplFig15 = ((synLTR_keep.sfs_group_mexnull_p + ggtitle('Filtered LTRs rooted by Mexicana')) | plot_spacer()) /
  ((synLTR_keep.sfs_group_folded_p + ggtitle('Filtered LTRs')) | (NAM.SNP.mexi_folded_p + ggtitle('Filtered SNPs'))) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 12, face = 'bold'))
pdf("supplFig15.pdf", width=9,height=6,pointsize=12, paper='special')
supplFig15
dev.off()


## Figure 1
layout <- "
#AA
BBB"
Fig1 =  pan_TE_bs_plot_pcnt + topSize_plot_SD +
  plot_layout(design = layout) + plot_layout(ncol = 1, heights = c(2, 3))
# plot out
pdf("Figure 1.pdf", width=10,height=8,pointsize=12, paper='special')
Fig1
dev.off()

## Figure 2
layout <- "
AABBBB
CCCDDD"
Fig2 = PCA_biplot_group + topVar_nor_p + trop_temp_diff_rank_plot4 + 
  trop_temp_size_cmb_plot + plot_layout(design = layout) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))
# plot out
pdf("Figure 2.pdf", width=10,height=8,pointsize=12, paper='special')
Fig2
dev.off()

## Figure 3
# draw parts and combine
pdf("Figure 3_part1.pdf", width=8,height=2.5,pointsize=12, paper='special')
act2 + phd3 + ina2
dev.off()

pdf("Figure 3_part2.pdf", width=6,height=8,pointsize=12, paper='special')
LTR_fam_SI_size_pval_tropbig_plot / LTR_fam_SI_size_pval_tempbig_plot
dev.off()

pdf("Figure 3_part3.pdf", width=6,height=5,pointsize=12, paper='special')
sizelist_LTR_ridge_plot_group2
dev.off()

pdf("Figure 3_part4.pdf", width=6,height=3,pointsize=12, paper='special')
LTR_status_diff_plot2
dev.off()

## Figure 4
Fig4 = (LTR_um_group_count_plot2 + (LTR_tropbig_tissue_plot2 + guides(fill = 'none')) + 
  plot_layout(guides = 'collect')) - Fig5_part23 + plot_layout(design = 'A\nB', heights = c(3, 6.5)) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))
# plot out
pdf("Figure 4.pdf", width=10,height=9.5,pointsize=12, paper='special')
Fig4
dev.off()


## Save suppl. files
write.table(synLTR %>% select(!c('genetdist', 'rrate')), file = "synLTR.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(synLTR %>% select(!c('genetdist', 'rrate')), file = "synLTR.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(solo, file = "soloLTR.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(LTR %>% select(!c('supfam', 'SO', 'sort_class', 'CHG_UM', 'genetdist', 
                              'rrate', 'rrate.f', 'aclen', 'ACR', 'AC_dist')), 
            file = "intactLTR.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(Fam_count %>% select(!c('sort_class', 'draw_size', 'non_um_count_temp', 
                              'non_um_count_trop', 'UMR_group_fisher_pval')), 
            file = "panTE_Fam_Stats.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
