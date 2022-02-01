## PopTEvo R code for expression analyses
## Shujun Ou (shujun.ou.1@gmail.com)
## 02/01/2022

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


library('dplyr')
if(!require(DESeq2)){BiocManager::install("DESeq2")}
if(!require(BiocParallel)){BiocManager::install("BiocParallel")}
library('ggplot2')
if(!require(coin)){install.packages("coin")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(multcompView)){install.packages("multcompView")}
register(MulticoreParam(2))
#library("pheatmap")
library("ComplexHeatmap")
library("RColorBrewer")
library(grid)
library("stringr")
library(patchwork)
library(scales)
library(EnvStats)
library("ggforce")
library("concaveman")
library(gridExtra)
library("ggbeeswarm")
library(tidyr)
library('ggsignif')

datadir="/Users/oushujun/Google\ Drive/My\ Drive/study/Maize\ research/NAM/TE/expression"
setwd(datadir)
load("./TE_fam_exp.RData")

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
cond_colors = c("Low remove Low diff"='#9BD9F9', "Low remove High diff"='#9DDEC2', "High remove Low diff"='#D2D296',
                "High remove High diff"='#FABEBB')

############################################
#### read data and prepare dataset ####
############################################
# read raw count data
TE_fam_exp_ori = read.table('family_sum_combined_counts_allNAM_15Jun21.txt', sep = '\t', header = T)
names(TE_fam_exp_ori) = sub('[A-Za-z0-9]+_X?', '', names(TE_fam_exp_ori))
names(TE_fam_exp_ori) = names(TE_fam_exp_ori) %>% sub("Oh7b", "Oh7B", .) %>% 
  sub("IL14H", "Il14H", .) %>% sub("MS71", "Ms71", .)
noquote(names(TE_fam_exp_ori))

load('../Fam_count.RData')
load('../LTR_fam_info.RData')
LTR_fam_info$condition = Fam_count$condition[match(LTR_fam_info$id, Fam_count$TE)]

# read library info
info = read.table('column_info.txt', sep = '\t', header = T)
info$name = sub('[A-Za-z0-9]+_X?', '', info$name)
info$name = info$name %>% sub("Oh7b", "Oh7B", .) %>% 
  sub("IL14H", "Il14H", .) %>% sub("MS71", "Ms71", .)
info$genome = info$genome %>% sub("Oh7b", "Oh7B", .) %>% 
  sub("IL14H", "Il14H", .) %>% sub("MS71", "Ms71", .)
info$group = as.factor(NAM$group[match(info$genome, NAM$genome)])
info$id = paste(info$tissue, info$group, sep = '_')
head(info)
unique(info$tissue)

## specify size factor using all mapped reads
lib_size = read.table('RNA_reads_mapped.txt', header= F, sep = '\t')
names(lib_size) = c('name', 'mapped_size')
lib_size$name = sub("Oh7b", "Oh7B", lib_size$name) %>% sub("IL14H", "Il14H", .) %>% 
  sub("MS71", "Ms71", .) %>% sub("Mo18*", "Mo18W", .)
lib_size$genome = sub('.*_', '', lib_size$name)
lib_size = lib_size[which(lib_size$genome != 'B73Ab10' & lib_size$genome != 'B73_AB10'), ] # remove Ab10
te_gene_sum = data.frame(name = names(colSums(TE_fam_exp_ori)), te_gene = colSums(TE_fam_exp_ori))
te_gene_sum$gene = t(TE_fam_exp_ori[rownames(TE_fam_exp_ori)=='Gene',])
te_gene_sum$name = sub("Oh7b", "Oh7B", te_gene_sum$name) %>% sub("IL14H", "Il14H", .) %>% 
  sub("MS71", "Ms71", .) %>% sub("Mo18*", "Mo18W", .)
lib_size$te_gene = te_gene_sum$te_gene[match(lib_size$name, te_gene_sum$name)]
lib_size$gene = te_gene_sum$gene[match(lib_size$name, te_gene_sum$name)]
lib_size$te = lib_size$te_gene - lib_size$gene
lib_size$te_ratio = lib_size$te/lib_size$te_gene
lib_size$tegenes_mappedratio = lib_size$te_gene/lib_size$mapped_size
lib_size$te_mappedratio = lib_size$te/lib_size$mapped_size
lib_size = lib_size %>% mutate(outliers = case_when(mapped_size >= 7.5e7 | mapped_size <= 2.5e7 ~ 'yes', TRUE ~ 'no'))

# check libraries
mean(lib_size$te_ratio)
mean(lib_size$tegenes_mappedratio)
mean(lib_size$te_mappedratio)
TE_fam_dds_sizeFactor = estimateSizeFactors(TE_fam_dds)$sizeFactor 
plot(lib_size$mapped_size[match(colnames(TE_fam_dds), lib_size$name)], TE_fam_dds_sizeFactor)
plot(lib_size$mapped_size, lib_size$te)
t.test(lib_size[lib_size$genome %in% NAM_trop,]$te_gene, lib_size[lib_size$genome %in% NAM_temp,]$te_gene)
lib_size = arrange(transform(lib_size, genome=factor(genome, levels=NAM_ID)), genome)
ggplot(lib_size, aes(genome, te_gene, fill = genome)) + geom_boxplot() + 
  scale_fill_manual(values = NAM_colors) +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95))
ggplot(lib_size, aes(te_gene, mapped_size, color = outliers)) + geom_point()
ggplot(lib_size, aes(mapped_size, te_ratio, color = outliers)) + geom_point()
tegenes_mappedratio_plot = ggplot(lib_size, aes(mapped_size, tegenes_mappedratio, color = outliers)) + geom_point()
te_ratio_plot = ggplot(lib_size, aes(te_gene, te_ratio, color = outliers)) + geom_point()
te_plot = ggplot(lib_size, aes(te, mapped_size, color = outliers)) + geom_point()
gene_plot = ggplot(lib_size, aes(gene, mapped_size, color = outliers)) + geom_point()

tegenes_mappedratio_plot + te_ratio_plot + te_plot + plot_annotation(tag_levels = 'A')

cor.test(lib_size$te_gene, lib_size$te_ratio)
cor.test(lib_size$te, lib_size$gene)
cor.test(lib_size$te, lib_size$mapped_size)
cor.test(lib_size$gene, lib_size$mapped_size)
cor.test(lib_size$te_gene, lib_size$mapped_size)
subset(lib_size, tegenes_mappedratio>0.5) #1 pair-end reads = 1 fragment in the raw count table

# keep only trop and temp
TE_fam_exp = TE_fam_exp_ori[, str_subset(names(TE_fam_exp_ori), 'M37W|Mo18W|Tx303', negate = T)]
info = subset(info, !genome %in% NAM_adx)

# keep pan TE fam only
TE_fam_exp = subset(TE_fam_exp, str_detect(rownames(TE_fam_exp), ':|\\(', negate = T))
TE_fam_exp[rownames(TE_fam_exp)=='Gene',] #check gene is still here
TE_fam_exp_ori[rownames(TE_fam_exp_ori)=='Gene',]

##################################################
### trop temp fam expressions for each tissue #######
##################################################
# keep only trop and temp
TE_fam_each_exp = TE_fam_exp_ori[, str_subset(names(TE_fam_exp_ori), 'M37W|Mo18W|Tx303', negate = T)]
TE_fam_each_exp = subset(TE_fam_each_exp, str_detect(rownames(TE_fam_each_exp), ':|\\(', negate = T))

# distinguish groups and tissues, construct deseq matrix
unique(info$id)
TE_fam_each_dds = DESeqDataSetFromMatrix(countData = TE_fam_each_exp,
                                    colData = subset(info, name %in% colnames(TE_fam_each_exp)),
                                    design = ~ id)

# add size factors
sizeFactors(TE_fam_each_dds) = sizeFactor_with_Genes4
mcols(TE_fam_each_dds)$basepairs = Fam_count$Size[match(rownames(TE_fam_each_exp), Fam_count$TE)]*1000000 # add TE fam length (all elements)

# run DESeq to get differentially expressed fams
TE_fam_each_dds <- DESeq(TE_fam_each_dds) #slow step

# get diff exp fams for each tissues and MA plots for each tissues
TE_fam_each_dds$id = factor(TE_fam_each_dds$id)
TE_fam_each_dds_sum <- vector("list", 1)
tissues = unique(info$tissue)

## expression level
# merge dataframes
TE_fam_each_dds_all <- Reduce(rbind, TE_fam_each_dds_sum)
TE_fam_each_dds_all$id = rownames(TE_fam_each_dds_all)
TE_fam_each_dds_all$tissue = factor(rep(tissues, each = length(unique(TE_fam_each_dds_all$id))))
TE_fam_each_dds_LTR = subset(TE_fam_each_dds_all, id %in% subset(Fam_count, str_detect(Supfam, 'LTR/'))$TE)
TE_fam_each_dds_LTR$condition = factor(Fam_count$condition[match(TE_fam_each_dds_LTR$id, Fam_count$TE)])
TE_fam_each_dds_LTR$status_class = factor(Fam_count$status_class[match(TE_fam_each_dds_LTR$id, Fam_count$TE)])
TE_fam_each_dds_LTR$Supfam = factor(Fam_count$Supfam[match(TE_fam_each_dds_LTR$id, Fam_count$TE)])

TE_fam_each_dds_all_sig = unique(rownames(subset(TE_fam_each_dds_all, padj<.05)))
TE_fam_each_dds_all_tropbig_sig = unique(rownames(subset(TE_fam_each_dds_all, padj<.05 & log2FoldChange>0)))
TE_fam_each_dds_all_tempbig_sig = unique(rownames(subset(TE_fam_each_dds_all, padj<.05 & log2FoldChange<0)))
LTR_trop_big_only = subset(trop_big_only, trop_big_only %in% Fam_count$TE[str_detect(Fam_count$Supfam, 'LTR/')])
LTR_temp_big_only = subset(temp_big_only, temp_big_only %in% Fam_count$TE[str_detect(Fam_count$Supfam, 'LTR/')])

tissue_sig = as.data.frame(table(subset(TE_fam_each_dds_all, padj<.05 & log2FoldChange>0)$tissue))
names(tissue_sig) = c('tissue', 'tropbig')
plot(tissue_sig$tropbig)
chisq.test(tissue_sig$tropbig)
library(EnvStats)
rosnerTest(tissue_sig$tropbig, k = 3) # no outliers
library(outliers)
grubbs.test(tissue_sig$tropbig) # no outlier
grubbs.test(tissue_sig$tropbig, opposite = T) # no outlier
tissue_sig_plot = ggplot(tissue_sig) + geom_col(aes(tissue, tropbig)) +
  labs(y = 'DE TE families (#)', x = 'Tissue') + theme_bw() +
  theme(axis.text.x = element_text(angle=35, vjust=1, hjust=0.95))
tissue_sig_plot

sum(subset(Fam_count, str_detect(Supfam, 'LTR'))$trop_temp_diff)
# 40.86417
sum(subset(Fam_count, TE %in% TE_fam_each_dds_all_sig)$trop_temp_diff)
# 23.00877 (mean ratio), 30.69902 (test each tissue)
sum(subset(Fam_count, TE %in% TE_fam_each_dds_all_tropbig_sig)$trop_temp_diff)
# 5.114903 (mapped), 15.42531 (mean ratio), 21.60929 (each tissue)
sum(subset(Fam_count, TE %in% TE_fam_each_dds_all_tempbig_sig)$trop_temp_diff)
# 22.92109 (mapped), 7.583454 (mean ratio), 12.53269 (each tissue)
sum(subset(Fam_count, TE %in% TE_fam_each_dds_all_tropbig_sig & trop_temp_diff>0)$trop_temp_diff)
# 17.31096 (mean ratio), 23.87067 (each tissue)
sum(subset(Fam_count, TE %in% TE_fam_each_dds_all_tempbig_sig & trop_temp_diff<0)$trop_temp_diff)
# -2.598052 (each tissue)
sum(subset(Fam_count, TE %in% TE_fam_each_dds_all_tropbig_sig & trop_temp_diff<0)$trop_temp_diff)
# -2.261381 (each tissue)
sum(subset(Fam_count, TE %in% TE_fam_each_dds_all_tempbig_sig & trop_temp_diff>0)$trop_temp_diff)
# 10.01374 (mean ratio), 15.13074 (each tissue)
sum(subset(Fam_count, TE %in% TE_fam_each_dds_all_tropbig_sig & category == 'Tropical amplification')$trop_temp_diff)
# 10.2569 (mean ratio), 16.08286 (each tissue)
sum(subset(Fam_count, TE %in% TE_fam_each_dds_all_tempbig_sig & category == 'Temperate removal')$trop_temp_diff)
# 4.426213 (mean ratio), 5.201592 (each tissue)

# find fams only sig exp in trop or temp
trop_big = unique(subset(TE_fam_each_dds_all, padj<.05 & log2FoldChange>0)$id) #1948
temp_big = unique(subset(TE_fam_each_dds_all, padj<.05 & log2FoldChange<0)$id) #1620
trop_big_only = trop_big[!(trop_big %in% temp_big)] #1789
temp_big_only = temp_big[!(temp_big %in% trop_big)] #1461

# find fams not conversely sig exp in tassel and anther
levels(TE_fam_each_dds_all$tissue)
trop_repro_high = unique(subset(TE_fam_each_dds_all,tissue %in% c('anther', 'tassel') & padj<.05 & log2FoldChange>0)$id) #925
temp_repro_high = unique(subset(TE_fam_each_dds_all,tissue %in% c('anther', 'tassel') & padj<.05 & log2FoldChange<0)$id) #590
trop_repro_high_only = setdiff(trop_repro_high, temp_repro_high) #916
trop_big_repro_consistent = trop_big[!(trop_big %in% temp_repro_high)] #1917
temp_big_repro_consistent = temp_big[!(temp_big %in% trop_repro_high)] #1566

sum(subset(Fam_count, TE %in% trop_big_only)$trop_temp_diff)
# 18.16633 (each tissue)
sum(subset(Fam_count, TE %in% trop_big_only & category == 'Tropical amplification')$trop_temp_diff)
# 13.66431 (each tissue)
sum(subset(Fam_count, TE %in% trop_big_repro_consistent & category == 'Tropical amplification')$trop_temp_diff)
# 16.07084 (each tissue)
# 16.08286/40.86417 = 0.3935687, 16.07084/40.86417 = 0.3932746
sum(subset(Fam_count, TE %in% temp_big_only)$trop_temp_diff)
# 9.089727 (each tissue)
sum(subset(Fam_count, TE %in% temp_big_only & category == 'Temperate removal')$trop_temp_diff)
# 4.812099 (each tissue)
sum(subset(Fam_count, TE %in% temp_big_repro_consistent & category == 'Temperate removal')$trop_temp_diff)
# 4.812099 (each tissue)

# focus on large contributon (>0.1Mb) families
trop_amp_large = as.character(subset(Fam_count, category == 'Tropical amplification' & trop_temp_diff > 0.1)$TE) # 50 fams
temp_rmv_large = as.character(subset(Fam_count, category == 'Temperate removal' & trop_temp_diff > 0.1)$TE) # 11 fams
diff0.1 = as.character(subset(Fam_count, trop_temp_diff > 0.1)$TE)
diff0.025 = as.character(subset(Fam_count, abs(trop_temp_diff) > 0.025)$TE)
diff0.025_LTR = as.character(subset(Fam_count, abs(trop_temp_diff) > 0.025 & 
                                      str_detect(Supfam, 'LTR/'))$TE)
length(subset(Fam_count, trop_temp_diff > 0.025)$TE)
length(subset(Fam_count, trop_temp_diff > 0.025)$TE)
trop_large0.025_LTR = as.character(subset(Fam_count, str_detect(Supfam, 'LTR/') & trop_temp_diff > 0.025)$TE) #138
temp_large0.025_LTR = as.character(subset(Fam_count, str_detect(Supfam, 'LTR/') & trop_temp_diff < -0.025)$TE) #27
trop_amp_large0.025 = as.character(subset(Fam_count, category == 'Tropical amplification' & trop_temp_diff > 0.025)$TE)
temp_rmv_large0.025 = as.character(subset(Fam_count, category == 'Temperate removal' & trop_temp_diff > 0.025)$TE)

trop_amp_large = trop_amp_large[order(match(trop_amp_large, Fam_count[order(Fam_count$trop_temp_diff, decreasing = T),]$TE))] # order from large diff to small diff
temp_rmv_large = temp_rmv_large[order(match(temp_rmv_large, Fam_count[order(Fam_count$trop_temp_diff, decreasing = T),]$TE))]
sum(subset(Fam_count, TE %in% trop_amp_large)$trop_temp_diff)
sum(subset(Fam_count, TE %in% temp_rmv_large)$trop_temp_diff)
sum(subset(Fam_count, str_detect(Supfam, 'LTR/'))$trop_temp_diff) #40.84318
sum(subset(Fam_count, TE %in% LTR_trop_big_only)$trop_temp_diff) #16.80755
16.80755/40.84318
sum(subset(Fam_count, TE %in% LTR_trop_big_only & category == 'Tropical amplification')$trop_temp_diff) #13.66431
sum(subset(Fam_count, TE %in% LTR_temp_big_only)$trop_temp_diff) #9.415073

sum(trop_amp_large %in% TE_fam_each_dds_all_tropbig_sig) #26 fams are significant higher exp in trop
sum(trop_amp_large %in% TE_fam_each_dds_all_tempbig_sig) #15 fams are significant higher exp in temp
sum(trop_amp_large0.025 %in% TE_fam_each_dds_all_tropbig_sig) #30 fams are significant higher exp in trop
sum(trop_amp_large0.025 %in% TE_fam_each_dds_all_tempbig_sig) #20 fams are significant higher exp in temp
fisher.test(rbind(c(30, 140), c(20, 150))) #% for fam diff >0.025,  trop amp trop higher vs trop amp temp higher, p-value = 0.1677
sum(trop_amp_large0.025 %in% LTR_trop_big_only) #23 LTR fams are significant higher exp in trop
sum(trop_amp_large0.025 %in% LTR_temp_big_only) #13 fams are significant higher exp in temp
sum(trop_large0.025_LTR %in% LTR_trop_big_only) #47 LTR trop 0.025larger fams are consistently and significantly higher exp in trop
sum(temp_large0.025_LTR %in% LTR_temp_big_only) #9 LTR temp 0.025larger fams are consistently and significantly higher exp in trop

# get lists
trop_large0.025 = subset(Fam_count, trop_temp_diff > 0.025)$TE #170
trop_large0.025_LTR #138
trop_amp = subset(Fam_count, category == 'Tropical amplification')$TE #131
temp_rmv = subset(Fam_count, category == 'Temperate removal')$TE #36
LTR_trop_big = subset(Fam_count, str_detect(Supfam, 'LTR/') & TE %in% trop_big)$TE #926
LTR_trop_big_only #839
trop_repro_high_only #916

# slice data
sum(subset(Fam_count, TE %in% intersect(intersect(trop_large0.025_LTR, LTR_trop_big_only), trop_amp))$trop_temp_diff) #13.48747
sum(subset(Fam_count, TE %in% trop_large0.025)$trop_temp_diff) #42.95978
sum(Fam_count$trop_temp_diff) #42.39669
fisher.test(rbind(c(23, 170-23), c(13.48747, 42.95978 - 13.48747))) #p-value = 0.01112

length(intersect(trop_large0.025_LTR, LTR_trop_big_only)) #47
length(intersect(intersect(trop_large0.025_LTR, LTR_trop_big_only), trop_amp)) #23
length(setdiff(trop_large0.025_LTR, LTR_trop_big_only)) #91
length(intersect(setdiff(trop_large0.025_LTR, LTR_trop_big_only), trop_amp)) #40
fisher.test(rbind(c(23, 24), c(40, 51))) #LTR_trop_big_only or not? trop_amp or not? not significant p-value = 0.5935

length(intersect(trop_large0.025_LTR, LTR_trop_big)) #62
length(intersect(intersect(trop_large0.025_LTR, LTR_trop_big), trop_amp)) #30
length(setdiff(trop_large0.025_LTR, LTR_trop_big)) #76
length(intersect(setdiff(trop_large0.025_LTR, LTR_trop_big), trop_amp)) #33
fisher.test(rbind(c(30, 32), c(33, 43))) #LTR_trop_big or not? trop_amp or not? not significant p-value = 0.5935

length(intersect(trop_amp, LTR_trop_big_only)) #45
length(intersect(intersect(trop_amp, LTR_trop_big_only), trop_large0.025_LTR)) #23
length(setdiff(trop_amp, LTR_trop_big_only)) #86
length(intersect(setdiff(trop_amp, LTR_trop_big_only), trop_large0.025_LTR)) #40
fisher.test(rbind(c(23, 22), c(40, 46))) #LTR_trop_big_only or not? diff0.025 or not? p-value = 0.7133

length(intersect(trop_large0.025_LTR, trop_repro_high_only)) #37
length(intersect(intersect(trop_large0.025_LTR, trop_repro_high_only), trop_amp)) #19
length(setdiff(trop_large0.025_LTR, trop_repro_high_only)) #101
length(intersect(setdiff(trop_large0.025_LTR, trop_repro_high_only), trop_amp)) #44
fisher.test(rbind(c(19, 18), c(44, 57))) #trop_repro_high_only or not? trop_amp or not? not significant p-value = 0.4454

sum(subset(Fam_count, TE %in% intersect(trop_large0.025_LTR, LTR_trop_big_only))$trop_temp_diff) #16.68064
sum(subset(Fam_count, TE %in% intersect(intersect(trop_large0.025_LTR, LTR_trop_big_only), trop_amp))$trop_temp_diff) #13.48747
sum(subset(Fam_count, TE %in% setdiff(trop_large0.025_LTR, LTR_trop_big_only))$trop_temp_diff) #24.42397
sum(subset(Fam_count, TE %in% intersect(setdiff(trop_large0.025_LTR, LTR_trop_big_only), trop_amp))$trop_temp_diff) #14.07601
fisher.test(rbind(c(13.48747, 3.19317), c(14.07601, 10.34796))) #LTR_trop_big_only or not? trop_amp or not? not significant p-value = 0.177

sum(subset(Fam_count, TE %in% intersect(trop_large0.025_LTR, LTR_trop_big))$trop_temp_diff) #19.93304
sum(subset(Fam_count, TE %in% intersect(intersect(trop_large0.025_LTR, LTR_trop_big), trop_amp))$trop_temp_diff) # 15.88297
sum(subset(Fam_count, TE %in% setdiff(trop_large0.025_LTR, LTR_trop_big))$trop_temp_diff) #21.17157
sum(subset(Fam_count, TE %in% intersect(setdiff(trop_large0.025_LTR, LTR_trop_big), trop_amp))$trop_temp_diff) #11.68051
fisher.test(rbind(c(15.88297, 4.05007), c(11.68051, 9.49106))) #LTR_trop_big or not? trop_amp or not? not significant p-value = 0.177

sum(subset(Fam_count, TE %in% intersect(trop_amp, LTR_trop_big_only))$trop_temp_diff) #13.66431
sum(subset(Fam_count, TE %in% intersect(intersect(trop_amp, LTR_trop_big_only), trop_large0.025_LTR))$trop_temp_diff) # 13.48747
sum(subset(Fam_count, TE %in% setdiff(trop_amp, LTR_trop_big_only))$trop_temp_diff) #14.28509
sum(subset(Fam_count, TE %in% intersect(setdiff(trop_amp, LTR_trop_big_only), trop_large0.025_LTR))$trop_temp_diff) #14.07601
fisher.test(rbind(c(13.48747, 0.17684), c(14.07601, 0.20908))) #LTR_trop_big_only or not? diff0.025 or not? not significant p-value = 1

sum(subset(Fam_count, TE %in% intersect(trop_large0.025_LTR, trop_repro_high_only))$trop_temp_diff) #13.01023
sum(subset(Fam_count, TE %in% intersect(intersect(trop_large0.025_LTR, trop_repro_high_only), trop_amp))$trop_temp_diff) # 11.3311
sum(subset(Fam_count, TE %in% setdiff(trop_large0.025_LTR, trop_repro_high_only))$trop_temp_diff) #28.09438
sum(subset(Fam_count, TE %in% intersect(setdiff(trop_large0.025_LTR, trop_repro_high_only), trop_amp))$trop_temp_diff) #16.23239
fisher.test(rbind(c(11.3311, 1.67913), c(16.23239, 11.86199))) #trop_repro_high_only or not? trop_amp or not? not significant p-value = 0.1559


intersect(trop_amp_large0.025, LTR_trop_big_only) #23
sum(subset(Fam_count, TE %in% intersect(trop_amp_large0.025, LTR_trop_big_only))$trop_temp_diff) #13.48747
13.48747/40.84318 #= 33.02258%
subset(Fam_count, TE %in% intersect(trop_amp_large0.025, LTR_trop_big_only))
sort(subset(Fam_count, category == 'Tropical amplification')$trop_temp_diff, decreasing = T)

sum(trop_amp_large %in% trop_big_only) #20 fams are only significant higher exp in trop
sum(trop_amp_large %in% temp_big_only) #9 fams are only significant higher exp in temp
sum(diff0.1 %in% TE_fam_each_dds_all_tropbig_sig) #40 fams are significant higher exp in trop
fisher.test(rbind(c(26, 24), c(40, 44))) #% of trop amp trop higher vs % of 0.1diff fam trop higher, p-value = 0.7213
fisher.test(rbind(c(26, 24), c(15, 35))) #% of trop amp trop higher vs % of trop amp trop lower, p-value = 0.04143
fisher.test(rbind(c(26, 24), c(1948, 18756))) #% of trop amp trop higher vs % of all fam trop higher, p-value = 2.928e-14
hist(subset(Fam_count, category == 'Temperate removal')$trop_temp_diff)
subset(TE_fam_dds_sum, id %in% trop_amp_large)
dim(subset(TE_fam_dds_sum, id %in% trop_amp_large & log2FoldChange > 0)) #31, trop amp fams with higher exp in trop
dim(subset(TE_fam_dds_sum, id %in% temp_rmv_large & log2FoldChange > 0))
length(unique(rownames(TE_fam_each_dds_all)))
length(TE_fam_each_dds_all_tropbig_sig)/length(unique(rownames(TE_fam_each_dds_all)))
TE_fam_each_dds_all_sig

# permutation test
TE_fam_each_dds_all$Supfam = Fam_count$Supfam[match(TE_fam_each_dds_all$id, Fam_count$TE)]
TE_fam_each_dds_all$category = Fam_count$category[match(TE_fam_each_dds_all$id, Fam_count$TE)]
TE_fam_each_dds_all$trop_temp_diff = Fam_count$trop_temp_diff[match(TE_fam_each_dds_all$id, Fam_count$TE)]
tropamp_test = subset(Fam_count, str_detect(Supfam, 'LTR/')) %>% 
  select(TE, Supfam, Size, trop_temp_diff, category) %>% 
  mutate(tropbigsig = case_when(TE %in% TE_fam_each_dds_all_tropbig_sig ~ '1', TRUE ~ '0'),
          tropamp = case_when(TE %in% trop_amp_large ~ '1', TRUE ~ '0'))

tropamp_test$tropbigsig = factor(tropamp_test$tropbigsig, ordered=TRUE)
tropamp_test$tropamp = factor(tropamp_test$tropamp, ordered=TRUE)
independence_test(tropamp ~ tropbigsig, data = tropamp_test) #p-value < 2.2e-16
table(tropamp_test$tropbigsig)
table(tropamp_test$tropamp)
matrix(data=c(26, 50, 926, 15317), nrow = 2)
fisher.test(matrix(data=c(26, 50, 926, 15317), nrow = 2)) #p-value = 5.482e-14, there are more tropbigsig fams in tropamp than other families
926/15317

subset(Fam_count, trop_temp_diff>0.025)$TE
table(subset(tropamp_test, TE %in% subset(Fam_count, trop_temp_diff>0.025)$TE)$tropbigsig)
fisher.test(matrix(data=c(26, 50, 62, 76), nrow = 2)) #p-value = 0.7748, there are more tropbigsig fams in tropamp than other families

# set up summary function
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
library(plyr)
library(dplyr)
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


# plot all TE fam exp comparisons
plot_fam_exp = function(fam_list){
  fam_list = fam_list
  TE_fam_large_exp_plot = list()
  for (i in 1:length(fam_list)){
    TE_fam = fam_list[i]
    print(TE_fam)
    geneCounts <- plotCounts(TE_fam_each_dds, gene = TE_fam, intgroup = c("group","tissue"), returnData = TRUE)
    geneCounts$group = as.factor(geneCounts$group)
    geneCounts_sum = summarySE(geneCounts, measurevar = 'count', groupvars = c('group', 'tissue'))
    TE_fam_large_exp_plot[[i]] = ggplot(geneCounts_sum, aes(x = as.numeric(group), y = count, color = tissue, group = tissue)) +
        geom_jitter(data = geneCounts, aes(x = as.numeric(group), y = count, color = tissue, group = tissue), 
                    size = 1, width = 0.15, alpha = 0.5) + 
      geom_errorbar(aes(ymin=count-se, ymax=count+se), width=.1, position=position_dodge(0.1)) +
      geom_line(position=position_dodge(0.1)) +
      geom_point(position=position_dodge(0.1)) +
      scale_x_continuous(breaks = c(1, 2), labels = c('temp', 'trop')) +
      scale_y_log10() + labs(y = 'Normalized count', x = '', title = TE_fam) +
      theme(legend.position = "none", axis.text.x = element_text(size=5), axis.text.y = element_text(size=5),
            axis.title.y = element_text(size=5), plot.title = element_text(size=3),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "grey"))
  }
  return(TE_fam_large_exp_plot)
}

TE_fam = 'cinful_zeon_AC209373_10201'
subset(TE_fam_each_dds_all, id == 'cinful_zeon_AC209373_10201')
geneCounts <- plotCounts(TE_fam_each_dds, gene = TE_fam, intgroup = c("group","tissue"), returnData = TRUE)
sample = ggplot(geneCounts, aes(x = group, y = count, color = tissue, group = tissue)) +
  scale_y_log10() + labs(y = 'Normalized count', x = '', title = TE_fam) +
  geom_jitter(size = 1, width = 0.15, alpha = 0.5) + geom_line() + guides(color=guide_legend(ncol=10)) +
  theme(legend.position = "top", axis.text.x = element_text(size=5), axis.text.y = element_text(size=5),
        axis.title.y = element_text(size=5), plot.title = element_text(size=3),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"))
ggsave("exp_plot_sample.pdf", sample, width = 12, height =8)


## combine sig fams from different tissues
LTR_fam_tissue = subset(TE_fam_each_dds_all, padj<.05)
LTR_fam_tissue$trop_temp_diff = Fam_count$trop_temp_diff[match(LTR_fam_tissue$id, Fam_count$TE)]
LTR_fam_tissue$status_class = Fam_count$status_class[match(LTR_fam_tissue$id, Fam_count$TE)]
LTR_fam_tissue$condition = Fam_count$condition[match(LTR_fam_tissue$id, Fam_count$TE)]
LTR_fam_tissue$category = Fam_count$category[match(LTR_fam_tissue$id, Fam_count$TE)]
LTR_fam_tissue$Supfam = Fam_count$Supfam[match(LTR_fam_tissue$id, Fam_count$TE)]
LTR_fam_tissue = as.data.frame(LTR_fam_tissue)

sum(subset(Fam_count, str_detect(Supfam, 'LTR/'))$trop_temp_diff) #40.84318
sum(distinct(LTR_fam_tissue, id, .keep_all= TRUE)$trop_temp_diff)
sum(distinct(subset(LTR_fam_tissue, str_detect(Supfam, 'LTR/') & trop_temp_diff > 0), id, .keep_all= TRUE)$trop_temp_diff) #31.60446
31.60446/40.84318 #sig higher expressed LTR fams explained 77.4% LTR size diff
sum(distinct(subset(LTR_fam_tissue, category == 'Tropical amplification' & trop_temp_diff > 0), id, .keep_all= TRUE)$trop_temp_diff) #19.80673
19.80673/40.84318 #sig higher expressed trop amp LTR fams explained 48.5% LTR size diff

# plot
LTR_fam_tissue_count = as.data.frame(table(table(LTR_fam_tissue$id)))
LTR_fam_tissue_count_p = ggplot(LTR_fam_tissue_count, aes(Var1, Freq)) + 
  geom_bar(stat = 'identity') + theme_bw() + 
  labs(x = 'Number of tissues', y = 'DE TE families (#)')
LTR_fam_tissue_count_p


##################################################
######## trop temp fam expressions #######
##################################################
# all tissues pooled
TE_fam_dds = DESeqDataSetFromMatrix(countData = TE_fam_exp,
                                    colData = subset(info, name %in% colnames(TE_fam_exp)),
                                    design = ~ group + tissue)
# add size factors
sizeFactors(TE_fam_dds) = sizeFactor_with_Genes4
mcols(TE_fam_dds)$basepairs = Fam_count$Size[match(rownames(TE_fam_exp), Fam_count$TE)]*1000000 # add TE fam length (all elements)
TE_fam_dds <- DESeq(TE_fam_dds) #slow step

## size-factor normalized fpkm
TE_fam_fpkm = as.data.frame(fpkm(TE_fam_dds, robust = TRUE)) #get normalized FPKM based on basepair length using DESeq
summary(colSums(TE_fam_fpkm[,1:441], na.rm=T))
TE_fam_fpkm$id = rownames(TE_fam_fpkm)
TE_fam_fpkm$category = Fam_count$category[match(TE_fam_fpkm$id, Fam_count$TE)]

#log transformation and add stats
TE_fam_fpkm_log2 = log2(TE_fam_fpkm + 1)
TE_fam_fpkm_log2$max = apply(TE_fam_fpkm_log2[,1:441], 1, max)
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

# stats
tapply(TE_fam_fpkm_log2$nam_mean, TE_fam_fpkm_log2$condition, mean)
tapply(TE_fam_fpkm_log2$nam_mean, TE_fam_fpkm_log2$category, mean)
table(TE_fam_fpkm_log2$condition)[1:4]/table(Fam_count$condition)[1:4]
table(TE_fam_fpkm_log2$category)/table(Fam_count$category)
tapply(TE_fam_fpkm_log2$nam_mean, TE_fam_fpkm_log2$condition, sum)
tapply(subset(Fam_count, condition != 'other')$Size, subset(Fam_count, condition != 'other')$condition, sum)
tapply(subset(Fam_count, category != 'Other')$Size, subset(Fam_count, category != 'Other')$category, sum)
tapply(subset(Fam_count, category != 'Other')$Size, subset(Fam_count, category != 'Other')$category, median)

# plot global sum
TE_fam_fpkm_log2_exp_pcnt = data.frame(table(subset(TE_fam_fpkm_log2, max>1)$condition)/table(Fam_count$condition))
names(TE_fam_fpkm_log2_exp_pcnt) = c('condition', 'exp_freq')
TE_fam_fpkm_log2_exp_pcnt_plot = ggplot(TE_fam_fpkm_log2_exp_pcnt) + 
  geom_bar(aes(condition, exp_freq), stat = 'identity') + theme_bw() +
  scale_x_discrete(breaks=c("High remove High diff","High remove Low diff",
                            "Low remove High diff","Low remove Low diff", "Other"),
                   labels=c("Temp removal","Trop-Temp balanced",
                            "Trop amplification","Random drifting", "Others")) +
  labs(y = 'Expressed LTR families') + scale_y_continuous(labels = percent) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=10, angle=35, vjust=1, hjust=0.95))
TE_fam_fpkm_log2_exp_pcnt_plot

#boxplot for fpkm
TE_fam_fpkm_log2_long = gather(TE_fam_fpkm_log2[,c(1:441, 448:451)], library, fpkm, 
              '16DAP_embryo_MN03101_B97':'V18_tassel_MN01062_B73', factor_key=TRUE)
TE_fam_fpkm_long = gather(TE_fam_fpkm, library, fpkm, 
                               '16DAP_embryo_MN03101_B97':'V18_tassel_MN01062_B73', factor_key=TRUE)
TE_fam_fpkm_long$Supfam = Fam_count$Supfam[match(TE_fam_fpkm_long$id, Fam_count$TE)]

# plot raw fpkm
t.test(subset(TE_fam_fpkm_long, category=='Tropical amplification' & id %in% focused_list)$fpkm,
       subset(TE_fam_fpkm_long, category=='Temperate removal' & id %in% focused_list)$fpkm) #p-value = 3.886e-15
t.test(subset(Fam_count, category=='Tropical amplification' & TE %in% focused_list)$Size,
       subset(Fam_count, category=='Temperate removal' & TE %in% focused_list)$Size) #p-value = 0.6124

TE_fam_fpkm_accu_main_plot = ggplot(subset(TE_fam_fpkm_long, category %in% c('Tropical amplification', 'Temperate removal') & id %in% focused_list), 
                                         aes(category, fpkm)) + geom_bar(size = 0.01, stat = 'summary', fun = 'sum') + 
  labs(x='', y='Accumulated FPKM') +
  theme(axis.text.x = element_text(angle=15, vjust=1, hjust=0.95),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

TE_fam_fpkm_mean_main_plot = ggplot(subset(TE_fam_fpkm_long, id %in% focused_list & 
                                                  category %in% c('Temperate removal', 'Tropical amplification')),
                                         aes(category, fpkm)) + 
  geom_bar(position = "dodge", stat = "summary", fun = "mean") + 
  geom_point(stat = "summary", fun = "median", shape = 23, color = 'red', fill = 'red') +
  geom_signif(annotations="P = 3.9e-15", y_position = 0.0029, xmin=1, xmax=2, size = 0.5, 
              vjust = -0.2, tip_length = c(0.0002, 0.0029)) + 
  coord_cartesian(ylim = c(0, 0.0033)) +
  labs(x='', y='Mean FPKM') +
  theme(axis.text.x = element_text(angle=15, vjust=1, hjust=0.95),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
TE_fam_fpkm_mean_main_plot

TE_fam_size_mean_main_plot = ggplot(subset(Fam_count, TE %in% focused_list & 
                                             category %in% c('Temperate removal', 'Tropical amplification')),
                                    aes(category, Size)) + 
  geom_bar(position = "dodge", stat = "summary", fun = "sum") + 
  labs(x='', y='Accumulated size (Mb)') +
  theme(axis.text.x = element_text(angle=15, vjust=1, hjust=0.95),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
TE_fam_size_mean_main_plot

####################
### make figures ###
####################
## suppl Fig 10
layout <- "
AAA
BBB
"

supplFig10 = (tissue_sig_plot + LTR_fam_tissue_count_p) - 
  (TE_fam_fpkm_accu_main_plot + TE_fam_size_mean_main_plot + TE_fam_fpkm_mean_main_plot) + 
  plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = 'bold'))

# plot out
pdf("suppl Fig10.pdf", width=8,height=7,pointsize=12, paper='special')
supplFig10
dev.off()

