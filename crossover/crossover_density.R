#!/usr/bin/env Rscript

# read recombination crossover info
NAM.rmap = read.table('./NAM.crossover.B73v5.pos', header=F)
colnames(NAM.rmap) = c('panel', 'chr', 'pos', 'marker')
#plot(density(NAM.rmap$pos, n=10000, na.rm=T))
NAM.rmap$panel = NAM.rmap$panel %>% sub("Oh7b", "Oh7B", .) %>% sub("Hp301", "HP301", .) %>% 
  sub("IL14H", "Il14H", .) %>% sub("MS71", "Ms71", .) %>% sub("Mo18*", "Mo18W", .)
NAM.rmap$chr = sub("chr", "", NAM.rmap$chr)
table(NAM.rmap$panel, NAM.rmap$chr)

# convert to recombination density
NAM.rrate = data.frame()
for (i in unique(NAM.rmap$panel)){
  print(i)
  this_panel = subset(NAM.rmap, panel == i)
  for (j in 1:10){
    this_panel_chr = subset(this_panel, chr == j)
    this_panel_chr_density = density(this_panel_chr$pos, n=max(this_panel_chr$pos)/25000, na.rm=T)
    this_rrate = this_panel_chr_density$y
    this_rrate_normalized = (this_rrate-min(this_rrate))/(max(this_rrate)-min(this_rrate))
    this_panel_chr_density_df = data.frame(panel=rep(i, length(this_rrate)), 
                                           chr=rep(j, length(this_rrate)), 
                                           pos=round(this_panel_chr_density$x, digits = 0),
                                           rrate=this_rrate_normalized)
    NAM.rrate = rbind(NAM.rrate, subset(this_panel_chr_density_df, pos>0))
  }
}

# test data and save
NAM.rrate$rrate.scale = scale(NAM.rrate$rrate) #standardization to mean=0, var=1
NAM.rrate$genome = sub('B73x', '', NAM.rrate$panel)

ggplot(subset(NAM.rrate, chr == '1')) + geom_line(aes(pos, rrate, color = panel)) + facet_wrap(~chr, ncol = 1)
ggplot(subset(NAM.rrate, chr == '1')) + geom_bar(aes(panel, mean(rrate)), stat = 'identity')
plot(tapply(NAM.rrate$rrate, NAM.rrate$panel, mean, na.rm=T))
plot(tapply(NAM.rrate$rrate.scale, NAM.rrate$panel, mean, na.rm=T))
t.test(NAM.rrate[NAM.rrate$genome %in% NAM_trop,]$rrate,
       NAM.rrate[NAM.rrate$geno %in% NAM_nontrop,]$rrate) #p-value = 1.286e-06
rm(NAM.rmap)

write.table(NAM.rrate, file = "NAM.rrate_v5.txt", append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
rm(NAM.rrate)
