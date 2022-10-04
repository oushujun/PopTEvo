## PopTEvo R code for expression analyses
## Shujun Ou (shujun.ou.1@gmail.com)
## 09/03/2022

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(DESeq2)){BiocManager::install("DESeq2")}
library('dplyr')
library(tidyr)
library('ggplot2')
library("RColorBrewer")
library(patchwork)
library(gridExtra)
library("stringr")
library(plyr)

if(!require(BiocParallel)){BiocManager::install("BiocParallel")}
register(MulticoreParam(2))
if(!require(coin)){install.packages("coin")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(multcompView)){install.packages("multcompView")}
library("ComplexHeatmap")
library(grid)
library(scales)
library(EnvStats)
library("ggforce")
library("concaveman")
library("ggbeeswarm")
library('ggsignif')

datadir='/Users/oushujun/My Drive/study/Maize research/NAM/TE/expression'
setwd(datadir)
load("./TE_fam_exp_clean.RData")

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

#######################################
#### read data and prepare dataset ####
#######################################
## read raw count data
TE_fam_count_ori = read.table('./NAM_expression/NAM_merged_Counts_7Jun22.txt', sep = '\t', header = T)

# format names and check
names(TE_fam_count_ori) = sub('[A-Za-z0-9]+_X?', '', names(TE_fam_count_ori))
names(TE_fam_count_ori) = names(TE_fam_count_ori) %>% sub("Oh7b", "Oh7B", .) %>% 
  sub("IL14H", "Il14H", .) %>% sub("MS71", "Ms71", .)
noquote(names(TE_fam_count_ori))
rownames(TE_fam_count_ori) = TE_fam_count_ori$id
TE_fam_count_ori$id = NULL
colSums(TE_fam_count_ori)

# read fam info
Fam_count = read.table('../panTE_Fam_Stats.txt', header = T, sep = '\t')
LTR_fam_info = read.table('../panLTR_Fam_Info.txt', header = T, sep = '\t')
LTR_fam_info$condition = Fam_count$condition[match(LTR_fam_info$id, Fam_count$TE)]

# read library info
info_ori = read.table('column_info.txt', sep = '\t', header = T)
info_ori$name = sub('[A-Za-z0-9]+_X?', '', info_ori$name)
info_ori$name = info_ori$name %>% sub("Oh7b", "Oh7B", .) %>% 
  sub("IL14H", "Il14H", .) %>% sub("MS71", "Ms71", .)
info_ori$genome = info_ori$genome %>% sub("Oh7b", "Oh7B", .) %>% 
  sub("IL14H", "Il14H", .) %>% sub("MS71", "Ms71", .)
info_ori$group = as.factor(NAM$group[match(info_ori$genome, NAM$genome)])
info_ori$id = paste(info_ori$tissue, info_ori$group, sep = '_')
head(info_ori)
unique(info_ori$tissue)
dim(info_ori)

# keep only trop and temp
TE_fam_count = TE_fam_count_ori[, str_subset(names(TE_fam_count_ori), 'M37W|Mo18W|Tx303', negate = T)]
info = subset(info_ori, !genome %in% NAM_adx)
info$group = factor(info$group) # remove unused groups
info$genome = factor(info$genome)

# keep pan TE fam only
TE_fam_count = subset(TE_fam_count, str_detect(rownames(TE_fam_count), ':|\\(', negate = T))
TE_fam_count[rownames(TE_fam_count)=='Gene',] #check gene is still here
dim(TE_fam_count) #15957, 441


#######################################
#### Check mapping rate and filter ####
#######################################
## specify size factor using all mapped reads
lib_size = read.table('RNA_reads_mapped.txt', header= F, sep = '\t')
names(lib_size) = c('name', 'mapped_size')
lib_size$name = sub("Oh7b", "Oh7B", lib_size$name) %>% sub("IL14H", "Il14H", .) %>% 
  sub("MS71", "Ms71", .) %>% sub("Mo18.*", "Mo18W", .)
lib_size$genome = sub('.*_', '', lib_size$name)
lib_size = lib_size[which(lib_size$genome != 'B73Ab10' & lib_size$genome != 'B73_AB10'), ] # remove Ab10
te_gene_sum = data.frame(name = names(colSums(TE_fam_count_ori)), te_gene = colSums(TE_fam_count_ori))
te_gene_sum$gene = t(TE_fam_count_ori[rownames(TE_fam_count_ori)=='Gene',])
te_gene_sum$name = sub("Oh7b", "Oh7B", te_gene_sum$name) %>% sub("IL14H", "Il14H", .) %>% 
  sub("MS71", "Ms71", .) %>% sub("Mo18.*", "Mo18W", .)
lib_size$te_gene = te_gene_sum$te_gene[match(lib_size$name, te_gene_sum$name)]
lib_size$gene = te_gene_sum$gene[match(lib_size$name, te_gene_sum$name)]
lib_size$te = lib_size$te_gene - lib_size$gene
lib_size$te_ratio = lib_size$te/lib_size$te_gene
lib_size$tegenes_mappedratio = lib_size$te_gene/lib_size$mapped_size
lib_size$te_mappedratio = lib_size$te/lib_size$mapped_size

# flag possible outliers
lib_size = lib_size %>% mutate(outliers = case_when(mapped_size >= 7.5e7 | mapped_size <= 2.5e7 ~ 'yes', TRUE ~ 'no'))

# check libraries
mean(lib_size$te_ratio)
mean(lib_size$te_mappedratio)
plot(lib_size$mapped_size, lib_size$te)
t.test(lib_size[lib_size$genome %in% NAM_trop,]$te_gene, lib_size[lib_size$genome %in% NAM_temp,]$te_gene)
lib_size = arrange(transform(lib_size, genome=factor(genome, levels=NAM_ID)), genome)
ggplot(lib_size, aes(genome, te_gene, fill = genome)) + geom_boxplot() + 
  scale_fill_manual(values = NAM_colors) +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95))
te_gene_size_plot = ggplot(lib_size, aes(te_gene, mapped_size, color = outliers)) + geom_point()
ggplot(lib_size, aes(mapped_size, te_ratio, color = outliers)) + geom_point()
tegenes_mappedratio_plot = ggplot(lib_size, aes(mapped_size, tegenes_mappedratio, color = outliers)) + geom_point()
te_ratio_plot = ggplot(lib_size, aes(te_gene, te_ratio, color = outliers)) + geom_point()
te_plot = ggplot(lib_size, aes(te, mapped_size, color = outliers)) + geom_point()
gene_plot = ggplot(lib_size, aes(gene, mapped_size, color = outliers)) + geom_point()

# complie plots together
te_gene_size_plot + tegenes_mappedratio_plot + te_ratio_plot + te_plot + plot_annotation(tag_levels = 'A')

cor.test(lib_size$te_gene, lib_size$te_ratio)
cor.test(lib_size$te, lib_size$gene)
cor.test(lib_size$te, lib_size$mapped_size)
cor.test(lib_size$gene, lib_size$mapped_size)
cor.test(lib_size$te_gene, lib_size$mapped_size)
subset(lib_size, tegenes_mappedratio>0.5) #1 pair-end reads = 1 fragment in the raw count table
mean(lib_size$tegenes_mappedratio) # ~0.3959703*2 read classification rate
hist(lib_size$tegenes_mappedratio)
hist(lib_size$te_ratio)

# filter out abnormal libraries
subset(lib_size, te_ratio < 0.01)
subset(lib_size, mapped_size/te_gene > 4)
subset(lib_size, tegenes_mappedratio < 0.3)
noquote(subset(lib_size, te_ratio == 0)$name)
remove_list =subset(lib_size, tegenes_mappedratio < 0.3)$name
TE_fam_count = TE_fam_count[, !names(TE_fam_count) %in% remove_list]
info = subset(info, !name %in% remove_list)
dim(TE_fam_count)
dim(info)

TE_fam_count[rownames(TE_fam_count)=='Gene',] #check gene is still here

## plot total expression and fam size
fam_exp = rowSums(TE_fam_count_ori)
Fam_count$read_count = fam_exp[Fam_count$TE]
cor.test(Fam_count$Size, Fam_count$read_count) #all TE fams: 0.4360764, p-value < 2.2e-16
cor.test(subset(Fam_count, str_detect(Fam_count$Supfam, 'LTR/'))$Size, 
         subset(Fam_count, str_detect(Fam_count$Supfam, 'LTR/'))$read_count) #all LTR fams: 0.4486134, p-value < 2.2e-16
cor.test(subset(Fam_count, str_detect(Fam_count$Supfam, 'LTR/') & category != 'Other')$Size, 
         subset(Fam_count, str_detect(Fam_count$Supfam, 'LTR/') & category != 'Other')$read_count) #classified LTR fams: 0.4471437, p-value < 2.2e-16

fam_size_read_count_plot = ggplot(subset(Fam_count, str_detect(Fam_count$Supfam, 'LTR/') & category != 'Other')) + 
  geom_point(aes(Size*1000, read_count, size = abs(trop_temp_diff), color = category), shape = 20) +
  scale_color_manual(breaks = c("Temperate removal", "Tropical amplification", 
                                "Balanced (trop > temp)", "Drifting (trop > temp)", 
                                "Tropical removal", "Temperate amplification",
                                "Balanced (trop < temp)", "Drifting (trop < temp)",
                                "Other"),
                     values = c('#FABEBB', '#9DDEC2', 'dark grey', 'dark grey', 
                                'dark grey', 'dark grey', 'dark grey', 'dark grey', 'dark grey'),
                     name = 'Classifications') +
  labs(x = 'Family size (kb)', y = 'Total RNA read count', size = 'Trop - Temp diff (Mb)')  + theme_bw() +
  scale_y_continuous(trans = 'log10') + scale_x_continuous(trans = 'log10')
fam_size_read_count_plot


#############################
#### Check data with PCA ####
#############################
# Build DESeq data including Gene expressions, use null design to find abnormal libraries
dds = DESeqDataSetFromMatrix(countData = TE_fam_count,
                             colData = subset(info, name %in% colnames(TE_fam_count)),
                             design = ~ 1)

# rlog normalization with blinding design to find abnormal libraries
vst <- varianceStabilizingTransformation(dds, blind = TRUE)
dat = plotPCA(vst, intgroup="genome", returnData = T)
dat$genome = info$genome[match(dat$name, info$name)]
dat$group = info$group[match(dat$name, info$name)]
dat$tissue = info$tissue[match(dat$name, info$name)]
proportion<-sprintf("%.2f", attr(dat, "percentVar")*100)
ggplot(aes(PC1, PC2, color=genome, shape=group),data=dat) + geom_point(size=3) +
  xlab(paste0('PC1', " (", proportion[1], "%)")) + ylab(paste0('PC2', " (", proportion[2], "%)"))
ggplot(aes(PC1, PC2, color=tissue, shape=group),data=dat) + geom_point(size=3) +
  xlab(paste0('PC1', " (", proportion[1], "%)")) + ylab(paste0('PC2', " (", proportion[2], "%)"))


# Build DESeq data including Gene expressions, use null design to find abnormal libraries
dds.bio = DESeqDataSetFromMatrix(countData = TE_fam_count,
                                    colData = subset(info, name %in% colnames(TE_fam_count)),
                                    design = ~ tissue + group + tissue:group)

# rlog normalization with experimental design to reveal biological patterns
vst.bio <- varianceStabilizingTransformation(dds.bio, blind = FALSE)
dat.bio = plotPCA(vst.bio, intgroup="genome", returnData = T)
dat.bio$genome = info$genome[match(dat.bio$name, info$name)]
dat.bio$group = info$group[match(dat.bio$name, info$name)]
dat.bio$tissue = info$tissue[match(dat.bio$name, info$name)]
proportion<-sprintf("%.2f", attr(dat.bio, "percentVar")*100)

# Plot PCA
mytheme = theme_bw() + theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size=10),
                             axis.title.y = element_text(size=12), axis.text.y = element_text(size=10),
                             legend.text = element_text(size=10), legend.title = element_blank())

RNA_libraries_tissue_group_PCA = ggplot(aes(PC1, PC2, color=tissue, shape=group),data=dat.bio) + geom_point(size=3) +
  xlab(paste0('PC1', " (", proportion[1], "%)")) + ylab(paste0('PC2', " (", proportion[2], "%)")) +
  ggtitle('PCA using RNA data of all shared TE families') + mytheme
RNA_libraries_tissue_group_PCA

#read count correlations
vst.norm.counts <- assay(vst.bio)
distance.vst <- as.dist(1-cor(vst.norm.counts, method="pearson")) 
tiff("vstTransformed_pearsonCorrelation.tiff", width = 24, height = 6, units = 'in', res = 300)  
plot(hclust(distance.vst), label=colnames(vst.norm.counts), cex = 0.3,
     main="vst transformed read counts\ndistance: Pearson correlation") 
dev.off()


################################
#### Normalize library size ####
################################
# set up a function to calculate geometric mean without 0 cells
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# use total TE read count as size factor
TE_count_total = as.numeric(lib_size$te)
names(TE_count_total) = lib_size$name
TE_count_total = subset(TE_count_total, names(TE_count_total) %in% names(TE_fam_count))
sizefactor = TE_count_total/gm_mean(TE_count_total)

## estimate group_tissue effects, construct deseq matrix
unique(info$id)
TE_fam_dds = DESeqDataSetFromMatrix(countData = TE_fam_count,
                                    colData = subset(info, name %in% colnames(TE_fam_count)),
                                    design = ~ id)

# add size factors
sizeFactors(TE_fam_dds) = sizefactor
mcols(TE_fam_dds)$basepairs = Fam_count$Size[match(rownames(TE_fam_count), Fam_count$TE)]*1000000 # add TE fam length (all elements)

# run DESeq to get differentially expressed fams
TE_fam_des = DESeq(TE_fam_dds, test = "Wald") #slow step


## estimate group effects, construct deseq matrix
TE_fam_dds = DESeqDataSetFromMatrix(countData = TE_fam_count,
                                    colData = subset(info, name %in% colnames(TE_fam_count)),
                                    design = ~ group + tissue + group:tissue)

# add size factors
sizeFactors(TE_fam_dds) = sizefactor
mcols(TE_fam_dds)$basepairs = Fam_count$Size[match(rownames(TE_fam_count), Fam_count$TE)]*1000000 # add TE fam length (all elements)

# run DESeq to get differentially expressed fams between groups
TE_fam_des_g = DESeq(TE_fam_dds, test = "Wald") #slow step


###################################
#### trop temp fam expressions ####
###################################
## size-factor normalized fpkm
TE_fam_fpkm = as.data.frame(fpkm(TE_fam_des, robust = TRUE)) #get normalized FPKM based on basepair length using DESeq
TE_fam_fpkm = subset(TE_fam_fpkm, rownames(TE_fam_fpkm) %in% Fam_count$TE) #keep just TE fam data
summary(colSums(TE_fam_fpkm[,1:433], na.rm=T))

#log transformation and add stats
TE_fam_fpkm_log2 = log2(TE_fam_fpkm + 1)
TE_fam_fpkm_log2$max = apply(TE_fam_fpkm_log2[,1:433], 1, max)
TE_fam_fpkm_log2$nam_mean = rowMeans(TE_fam_fpkm_log2[, grep(pattern = paste(names(NAM_colors), collapse="|"), 
                                                             x = names(TE_fam_fpkm_log2))], na.rm = T)
TE_fam_fpkm_log2$trop_mean = rowMeans(TE_fam_fpkm_log2[, grep(pattern = paste(NAM_trop, collapse="|"), 
                                                              x = names(TE_fam_fpkm_log2))], na.rm = T)
TE_fam_fpkm_log2$temp_mean = rowMeans(TE_fam_fpkm_log2[, grep(pattern = paste(NAM_nontrop, collapse="|"), 
                                                              x = names(TE_fam_fpkm_log2))], na.rm = T)
TE_fam_fpkm_log2$trop_se = apply(TE_fam_fpkm_log2[, grep(pattern = paste(NAM_trop, collapse="|"), 
                                                         x = names(TE_fam_fpkm_log2))], 1, function(x) sd(x)/sqrt(length(x)))
TE_fam_fpkm_log2$temp_se = apply(TE_fam_fpkm_log2[, grep(pattern = paste(NAM_nontrop, collapse="|"), 
                                                         x = names(TE_fam_fpkm_log2))], 1, function(x) sd(x)/sqrt(length(x)))
TE_fam_fpkm_log2$id = rownames(TE_fam_fpkm_log2)
TE_fam_fpkm_log2$Supfam = Fam_count$Supfam[match(TE_fam_fpkm_log2$id, Fam_count$TE)]
TE_fam_fpkm_log2$condition = Fam_count$condition[match(TE_fam_fpkm_log2$id, Fam_count$TE)]
TE_fam_fpkm_log2$category = Fam_count$category[match(TE_fam_fpkm_log2$id, Fam_count$TE)]

# make long-format for fpkm data
TE_fam_fpkm_log2_long = gather(TE_fam_fpkm_log2[,c(1:433, 440:443)], library, fpkm, 
                               '16DAP_embryo_MN01101_B73':'V18_tassel_MN27062_Tzi8', factor_key=TRUE)
TE_fam_fpkm$id = rownames(TE_fam_fpkm)
TE_fam_fpkm$category = Fam_count$category[match(TE_fam_fpkm$id, Fam_count$TE)]
TE_fam_fpkm_long = gather(TE_fam_fpkm, library, fpkm, 
                          '16DAP_embryo_MN01101_B73':'V18_tassel_MN27062_Tzi8', factor_key=TRUE)
TE_fam_fpkm_long$Supfam = Fam_count$Supfam[match(TE_fam_fpkm_long$id, Fam_count$TE)]

# stat expression level
tapply(TE_fam_fpkm_long$fpkm, TE_fam_fpkm_long$category, sum) #1794.89163/86.10461 = 20.84548
tapply(Fam_count$Size, Fam_count$category, sum)
tapply(TE_fam_fpkm_long$fpkm, TE_fam_fpkm_long$category, sum)/ #normalize with total seq length
  tapply(Fam_count$Size, Fam_count$category, sum) #4.3038038/0.8877386 = 4.848053
t.test(subset(TE_fam_fpkm_long, category=='Tropical amplification')$fpkm,
       subset(TE_fam_fpkm_long, category=='Temperate removal')$fpkm) #p-value < 2.2e-16, trop amp larger


# plot raw fpkm
TE_fam_fpkm_accu_plot = ggplot(subset(TE_fam_fpkm_long, category %in% c('Tropical amplification', 'Temperate removal')), 
                               aes(category, fpkm)) + geom_bar(size = 0.01, stat = 'summary', fun = 'sum') + 
  labs(x='', y='Accumulated FPKM') +
  theme(axis.text.x = element_text(angle=15, vjust=1, hjust=0.95),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
TE_fam_fpkm_accu_plot

TE_fam_size_mean_plot = ggplot(subset(Fam_count, category %in% c('Temperate removal', 'Tropical amplification')),
                               aes(category, Size)) + 
  geom_bar(position = "dodge", stat = "summary", fun = "sum") + 
  labs(x='', y='Accumulated size (Mb)') +
  theme(axis.text.x = element_text(angle=15, vjust=1, hjust=0.95),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
TE_fam_size_mean_plot

TE_fam_fpkm_mean_plot = ggplot(subset(TE_fam_fpkm_long, category %in% c('Temperate removal', 'Tropical amplification')),
                               aes(category, fpkm)) + 
  geom_bar(position = "dodge", stat = "summary", fun = "mean") + 
  geom_point(stat = "summary", fun = "median", shape = 23, color = 'red', fill = 'red') +
  geom_signif(annotations="P < 1.0e-10", y_position = 0.032, xmin=1, xmax=2, size = 0.5, 
              vjust = -0.2, tip_length = c(0.0006, 0.00005)) + 
  coord_cartesian(ylim = c(0, 0.035)) +
  labs(x='', y='Mean FPKM') +
  theme(axis.text.x = element_text(angle=15, vjust=1, hjust=0.95),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
TE_fam_fpkm_mean_plot

TE_fam_fpkm_vio_plot = ggplot(subset(TE_fam_fpkm_long, 
                                     category %in% c('Temperate removal', 'Tropical amplification')),
                              aes(category, fpkm)) + 
  geom_violin(position = position_dodge(width = 1)) + 
  geom_point(position = position_dodge(width = 1), alpha = 0.5, width = 0.2) +
  labs(x='', y='Mean FPKM') +
  theme(axis.text.x = element_text(angle=15, vjust=1, hjust=0.95),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
TE_fam_fpkm_vio_plot


###################################################
#### trop temp fam expressions for each tissue ####
###################################################
# compare each tissue, log2FoldChange = log2(temp/trop), negatve value means expression higher in trop, p values are adjusted by Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
TE_fam_des_anther <- data.frame(TE = rownames(TE_fam_des), class = 'anther', results(TE_fam_des, contrast = c('id', 'anther_temp', 'anther_trop'), pAdjustMethod = "BH"))
TE_fam_des_ear <- data.frame(TE = rownames(TE_fam_des), class = 'ear', results(TE_fam_des, contrast = c('id', 'ear_temp', 'ear_trop'), pAdjustMethod = "BH"))
TE_fam_des_embryo <- data.frame(TE = rownames(TE_fam_des), class = 'embryp', results(TE_fam_des, contrast = c('id', 'embryo_temp', 'embryo_trop'), pAdjustMethod = "BH"))
TE_fam_des_endosperm <- data.frame(TE = rownames(TE_fam_des), class = 'endosperm', results(TE_fam_des, contrast = c('id', 'endosperm_temp', 'endosperm_trop'), pAdjustMethod = "BH"))
TE_fam_des_middle <- data.frame(TE = rownames(TE_fam_des), class = 'middle', results(TE_fam_des, contrast = c('id', 'middle_temp', 'middle_trop'), pAdjustMethod = "BH"))
TE_fam_des_root <- data.frame(TE = rownames(TE_fam_des), class = 'root', results(TE_fam_des, contrast = c('id', 'root_temp', 'root_trop'), pAdjustMethod = "BH"))
TE_fam_des_shoot <- data.frame(TE = rownames(TE_fam_des), class = 'shoot', results(TE_fam_des, contrast = c('id', 'shoot_temp', 'shoot_trop'), pAdjustMethod = "BH"))
TE_fam_des_tassel <- data.frame(TE = rownames(TE_fam_des), class = 'tassel', results(TE_fam_des, contrast = c('id', 'tassel_temp', 'tassel_trop'), pAdjustMethod = "BH"))
TE_fam_des_tip <- data.frame(TE = rownames(TE_fam_des), class = 'tip', results(TE_fam_des, contrast = c('id', 'tip_temp', 'tip_trop'), pAdjustMethod = "BH"))
TE_fam_des_base <- data.frame(TE = rownames(TE_fam_des), class = 'base', results(TE_fam_des, contrast = c('id', 'base_temp', 'base_trop'), pAdjustMethod = "BH"))

# compare groups, log2FoldChange = log2(temp/trop), negatve value means expression higher in trop
TE_fam_des_group = data.frame(TE = rownames(TE_fam_des), class = 'group', results(TE_fam_des_g, contrast = c('group', 'temp', 'trop'), pAdjustMethod = "BH"))

# combine results
TE_fam_res = rbind(TE_fam_des_anther, TE_fam_des_base, TE_fam_des_ear, TE_fam_des_embryo, 
                   TE_fam_des_endosperm, TE_fam_des_middle, TE_fam_des_root, TE_fam_des_shoot, 
                   TE_fam_des_tassel, TE_fam_des_tip, TE_fam_des_group)

# summarize results
sum(subset(Fam_count, TE %in% unique(subset(TE_fam_res, padj <= 0.05)$TE))$trop_temp_diff) #22.88565
sum(subset(Fam_count, TE %in% subset(TE_fam_res, padj <= 0.05 & class == 'group' & log2FoldChange < 0)$TE)$trop_temp_diff) #5.188679
sum(subset(Fam_count, TE %in% subset(TE_fam_res, padj <= 0.05 & log2FoldChange < 0)$TE)$trop_temp_diff) #19.39999

## expression level
TE_fam_res_LTR = subset(TE_fam_res, TE %in% subset(Fam_count, str_detect(Supfam, 'LTR/'))$TE)
TE_fam_res_LTR$condition = factor(Fam_count$condition[match(TE_fam_res_LTR$TE, Fam_count$TE)])
TE_fam_res_LTR$status_class = factor(Fam_count$status_class[match(TE_fam_res_LTR$TE, Fam_count$TE)])
TE_fam_res_LTR$Supfam = factor(Fam_count$Supfam[match(TE_fam_res_LTR$TE, Fam_count$TE)])

# find fams only sig exp in trop or temp
trop_big = unique(subset(TE_fam_res, padj<.05 & log2FoldChange<0)$TE) #3286
temp_big = unique(subset(TE_fam_res, padj<.05 & log2FoldChange>0)$TE) #1460
trop_big_only = trop_big[!(trop_big %in% temp_big)] #2927
temp_big_only = temp_big[!(temp_big %in% trop_big)] #1101
LTR_trop_big_only = subset(trop_big_only, trop_big_only %in% Fam_count$TE[str_detect(Fam_count$Supfam, 'LTR/')]) #1581
LTR_temp_big_only = subset(temp_big_only, temp_big_only %in% Fam_count$TE[str_detect(Fam_count$Supfam, 'LTR/')]) #634
length(LTR_trop_big_only) #1581

# find tissue outliers for fams higher exp in trop
tissue_sig = as.data.frame(table(subset(TE_fam_res, padj<.05 & log2FoldChange<0)$class))
names(tissue_sig) = c('tissue', 'tropbig')
plot(tissue_sig$tropbig)
chisq.test(tissue_sig$tropbig)
library(EnvStats)
rosnerTest(tissue_sig$tropbig, k = 3) # no outlier
library(outliers)
grubbs.test(tissue_sig$tropbig) # root is an outlier
grubbs.test(tissue_sig$tropbig, opposite = T) # no outlier
tissue_sig_plot = ggplot(subset(tissue_sig, tissue != 'group')) + geom_col(aes(tissue, tropbig)) +
  labs(y = 'DE TE families (#)', x = 'Tissue') + theme_bw() +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95))
tissue_sig_plot

sum(subset(Fam_count, str_detect(Supfam, 'LTR'))$trop_temp_diff)
# 34.42639
sum(subset(Fam_count, TE %in% LTR_trop_big_only & trop_temp_diff>0)$trop_temp_diff) #133 fams
# 20.49653/34.42639 = 0.5953726
sum(subset(Fam_count, TE %in% LTR_temp_big_only & trop_temp_diff<0)$trop_temp_diff)
# -1.332725/34.42639 = -0.03871231
sum(subset(Fam_count, TE %in% LTR_trop_big_only & trop_temp_diff>0 & category == 'Tropical amplification')$trop_temp_diff) #28 fams
# 11.01076/34.42639 = 0.3198349
sum(subset(Fam_count, TE %in% LTR_trop_big_only & trop_temp_diff>0 & category == 'Temperate removal')$trop_temp_diff) #9 fams
# 4.213802/34.42639 = 0.1224003
sum(subset(Fam_count, TE %in% LTR_temp_big_only & trop_temp_diff<0 & category == 'Temperate removal')$trop_temp_diff) #0 fams
# 0/34.42639 = 0
sum(subset(Fam_count, trop_temp_diff>0 & category == 'Tropical amplification' & 
             TE %in% LTR_trop_big_only & um_count > 0)$trop_temp_diff) # 9.961229, 19 fams
# 9.961229/34.42639 = 0.2893486


# randomly sample length(LTR_trop_big_only) families for 1000 times and see how much trop_temp_diff they contribute
trop_amp_rand_diff = 0
temp_rmv_rand_diff = 0
resample <- function(x, ...) x[sample.int(length(x), ...)]
for (i in 1:1000){
  random = resample(unique(TE_fam_res$TE), length(LTR_trop_big_only))
  trop_amp_rand = random[random %in% subset(Fam_count, trop_temp_diff>0 & category == 'Tropical amplification')$TE]
  temp_rmv_rand = random[random %in% subset(Fam_count, trop_temp_diff>0 & category == 'Temperate removal')$TE]
  trop_amp_rand_diff = trop_amp_rand_diff + sum(subset(Fam_count, TE %in% trop_amp_rand)$trop_temp_diff)
  temp_rmv_rand_diff = temp_rmv_rand_diff + sum(subset(Fam_count, TE %in% temp_rmv_rand)$trop_temp_diff)
}
trop_amp_rand_diff/1000 #1.948008
temp_rmv_rand_diff/1000 #1.052625

# construct a fisher test 
sum(subset(Fam_count, trop_temp_diff>0 & category == 'Tropical amplification' & 
             TE %in% LTR_trop_big_only)$trop_temp_diff) # 11.01076
sum(subset(Fam_count, trop_temp_diff>0 & category == 'Tropical amplification' & 
             TE %in% unique(TE_fam_res$TE))$trop_temp_diff) #19.21384
sum(subset(Fam_count, trop_temp_diff>0 & category == 'Temperate removal' & 
             TE %in% LTR_trop_big_only)$trop_temp_diff) # 4.213802
sum(subset(Fam_count, trop_temp_diff>0 & category == 'Temperate removal' & 
             TE %in% unique(TE_fam_res$TE))$trop_temp_diff) #10.96527
# null hypothesis: LTR_trop_big_only is a random list for trop amp families
fisher.test(rbind(c(1.948008, 11.01076), c(19.21384, 19.21384))) #p-value = 0.04841, reject the null hypothesis
11.01076/1.948008 #5.652318 times larger than random expectation
# null hypothesis: LTR_trop_big_only is a random list for temp rmv families
fisher.test(rbind(c(1.052625, 4.213802), c(10.96527, 10.96527))) #p-value = 0.3419, fail to reject the null hypothesis

#### plot expression norms (fam exp b/t trop and temp genomes)
## set up summary function
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
#install.packages("~/Downloads/doBy_4.6.11.tar.gz", repos = NULL, type="source")
# doBy package download: https://cran.r-project.org/web/packages/doBy/index.html
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

# a function to plot TE expression comparisons between trop and temp among all libraries
plot_fam_exp = function(fam_list){
  fam_list = fam_list
  TE_fam_large_exp_plot = list()
  for (i in 1:length(fam_list)){
    TE_fam = fam_list[i]
    print(TE_fam)
    geneCounts <- plotCounts(TE_fam_des, gene = TE_fam, intgroup = c("group","tissue"), returnData = TRUE)
    geneCounts$group = as.factor(geneCounts$group)
    geneCounts_sum = summarySE(geneCounts, measurevar = 'count', groupvars = c('group', 'tissue'))
    TE_fam_large_exp_plot[[i]] = ggplot(geneCounts_sum, aes(x = as.numeric(group), y = count, color = tissue, group = tissue)) +
      geom_jitter(data = geneCounts, aes(x = as.numeric(group), y = count, color = tissue, group = tissue), 
                  size = 1, width = 0.15, alpha = 0.5) + 
      geom_errorbar(aes(ymin=count-se, ymax=count+se), width=.1, position=position_dodge(0.1)) +
      geom_line(position=position_dodge(0.1)) +
      geom_point(position=position_dodge(0.1)) +
      scale_x_continuous(breaks = c(2, 3), labels = c('Temp', 'Trop')) +
      scale_y_log10() + labs(y = 'Normalized count', x = '', title = TE_fam) +
      theme(legend.position = "none", axis.text.x = element_text(size=5), axis.text.y = element_text(size=5),
            axis.title.y = element_text(size=5), plot.title = element_text(size=3),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "grey"))
  }
  return(TE_fam_large_exp_plot)
}

# plot an example to examine and also get the legend
TE_fam = 'cinful_zeon_AC209373_10201'
subset(TE_fam_res, TE == 'cinful_zeon_AC209373_10201')
geneCounts <- plotCounts(TE_fam_des, gene = TE_fam, intgroup = c("group","tissue"), returnData = TRUE)
sample = plot_fam_exp(TE_fam)[[1]] + theme(legend.position = "top") + guides(color=guide_legend(ncol=10))
ggsave("exp_plot_sample.pdf", sample, width = 12, height =8)


# identify large diff families and their expression significance levels
trop_amp_large = as.character(subset(Fam_count, category == 'Tropical amplification' & trop_temp_diff > 0.1)$TE) # 38 fams
temp_rmv_large = as.character(subset(Fam_count, category == 'Temperate removal' & trop_temp_diff > 0.1)$TE) # 14 fams
trop_amp_large = subset(Fam_count, TE %in% trop_amp_large)$TE[order(subset(Fam_count, TE %in% trop_amp_large)$trop_temp_diff, decreasing = T)] #order based on trop_temp_diff
temp_rmv_large = subset(Fam_count, TE %in% temp_rmv_large)$TE[order(subset(Fam_count, TE %in% temp_rmv_large)$trop_temp_diff, decreasing = T)] #order based on trop_temp_diff
trop_amp_large[trop_amp_large %in% LTR_trop_big_only] # *marked on Suppl. Fig. 11

# plot expression of trop amp and temp rmv fams
Trop_amp_large_exp_plot = plot_fam_exp(trop_amp_large)
Trop_amp_large_exp_plot_all = marrangeGrob(Trop_amp_large_exp_plot, ncol = 10, nrow = 5)
ggsave("Trop_amp_large_exp_plot_all_3.pdf", Trop_amp_large_exp_plot_all, width = 12, height =8)
temp_rmv_large_exp_plot = plot_fam_exp(temp_rmv_large)
temp_rmv_large_exp_plot_all = marrangeGrob(temp_rmv_large_exp_plot, ncol = 10, nrow = 2)
ggsave("temp_rmv_large_exp_plot_all2.pdf", temp_rmv_large_exp_plot_all, width = 12, height =3)


## plot tissue frequency of sig fams
LTR_fam_tissue_count = as.data.frame(table(table(subset(TE_fam_res, padj<.05 & class != 'group')$TE)))
LTR_fam_tissue_count_p = ggplot(LTR_fam_tissue_count, aes(Var1, Freq)) + 
  geom_bar(stat = 'identity') + theme_bw() + 
  labs(x = 'Number of tissues', y = 'DE TE families (#)')
LTR_fam_tissue_count_p


####################
### make figures ###
####################
## suppl Fig 10
layout <- "
AAA
BBB
"

supplFig10 = (tissue_sig_plot + LTR_fam_tissue_count_p) - 
  (TE_fam_fpkm_accu_plot + TE_fam_size_mean_plot + TE_fam_fpkm_mean_plot) + 
  plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = 'bold'))

# plot out
pdf("suppl Fig10.pdf", width=8,height=7,pointsize=12, paper='special')
supplFig10
dev.off()

# New suppl fig
pdf("suppl Fig new.pdf", width=8,height=6,pointsize=12, paper='special')
fam_size_read_count_plot
dev.off()

# save data
save(TE_fam_res_LTR, file = "LTR_fam_tissue.RData")

