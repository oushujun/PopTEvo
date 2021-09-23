#!/usr/bin/env Rscript

# Yinjie Qiu (qiuxx221@umn.edu)
# 05/24/2021

order = c("B73AB10", "B97", "CML103", "CML228", "CML247", "CML277", "CML322",
          "CML333", "CML52", "CML69", "HP301", "Il14H", "Ki11", "Ki3", "Ky21", "M162W",
          "M37W", "Mo18W", "Ms71", "NC350", "NC358", "Oh43", "Oh7b", "P39", "Tx303", "Tzi8")

for (i in order){
  B73 <- read.csv(paste0(i,"_B73_all_pair.txt"),sep="\t")
  B73AB <- read.csv(paste0(i,"_B73AB10_all_pair.txt"),sep="\t")
  B97 <- read.csv(paste0(i,"_B97_all_pair.txt"), sep = "\t")
  Ky21 <- read.csv(paste0(i,"_Ky21_all_pair.txt"), sep = "\t")
  M162W <- read.csv(paste0(i,"_M162W_all_pair.txt"), sep = "\t")
  Ms71 <- read.csv(paste0(i,"_Ms71_all_pair.txt"), sep = "\t")
  Oh7B <- read.csv(paste0(i,"_Oh7B_all_pair.txt"), sep = "\t")
  Oh43 <- read.csv(paste0(i,"_Oh43_all_pair.txt"), sep = "\t")
  M37W <- read.csv(paste0(i,"_M37W_all_pair.txt"), sep = "\t")
  Mo18W <- read.csv(paste0(i,"_Mo18W_all_pair.txt"), sep = "\t")
  Tx303 <- read.csv(paste0(i,"_Tx303_all_pair.txt"), sep = "\t")
  HP301 <- read.csv(paste0(i,"_HP301_all_pair.txt"), sep = "\t")
  Il14H <- read.csv(paste0(i,"_Il14H_all_pair.txt"), sep = "\t")
  P39 <- read.csv(paste0(i,"_P39_all_pair.txt"), sep = "\t")
  CML52 <- read.csv(paste0(i,"_CML52_all_pair.txt"), sep = "\t")
  CML69 <- read.csv(paste0(i,"_CML69_all_pair.txt"), sep = "\t")
  CML103 <- read.csv(paste0(i,"_CML103_all_pair.txt"), sep = "\t")
  CML228 <- read.csv(paste0(i,"_CML228_all_pair.txt"), sep = "\t")
  CML247 <- read.csv(paste0(i,"_CML247_all_pair.txt"), sep = "\t")
  CML277 <- read.csv(paste0(i,"_CML277_all_pair.txt"), sep = "\t")
  CML322 <- read.csv(paste0(i,"_CML322_all_pair.txt"), sep = "\t")
  CML333 <- read.csv(paste0(i,"_CML333_all_pair.txt"), sep = "\t")
  Ki3 <- read.csv(paste0(i,"_Ki3_all_pair.txt"), sep = "\t")
  Ki11 <- read.csv(paste0(i,"_Ki11_all_pair.txt"), sep = "\t")
  NC350 <- read.csv(paste0(i,"_NC350_all_pair.txt"), sep = "\t")
  NC358 <- read.csv(paste0(i,"_NC358_all_pair.txt"), sep = "\t")
  Tzi8 <- read.csv(paste0(i,"_Tzi8_all_pair.txt"), sep = "\t")
  anchor <- read.csv(paste0(i,"_anchor.txt"), sep = "\t")
  left_join(anchor, B73) %>% left_join(B73AB) %>% left_join(Tzi8) %>% 
    left_join(Ky21) %>% left_join(M162W) %>% left_join(Ms71) %>% left_join(Oh7B) %>% left_join(Oh43) %>% 
    left_join(M37W) %>% left_join(Mo18W) %>% left_join(NC350) %>% left_join(HP301) %>% left_join(Il14H) %>% 
    left_join(P39) %>% left_join(CML52) %>% left_join(CML69) %>% left_join(Ki11) %>% left_join(CML228) %>% 
    left_join(CML247) %>% left_join(CML277) %>% left_join(CML322) %>% left_join(CML333) %>% left_join(Ki3) %>%
    left_join(CML103) %>% left_join(Tx303) %>% left_join(NC358) %>% left_join(B97) %>% write.csv(file=paste0(i,"_pan_TE_pair.csv"))
}
