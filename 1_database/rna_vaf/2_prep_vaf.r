library(data.table)
library(tidyverse)
library(gridExtra)

vafs <- fread("/home/jusset/vaf_clean/data/vaf_ready.csv") %>% select(-V1)

saveRDS(vafs %>% filter(tier == "PANEL") %>% arrange(gene) %>% pull(gene) %>% unique(), 
        file = "/home/jusset/vaf_clean/data/panel_gene.Rds")

vaf_ready <- 
vafs %>% 
  filter( subclonal < .1, purple_cn > .5, !is.na(purple_af) ) %>% 
  mutate( purple_af_mx = ifelse(biallelic, 1, purple_af)) %>% #### changed to biallelic
  mutate( rna_vaf_est = ifelse(rna_vaf/purple_af_mx > 1, 1, rna_vaf/purple_af_mx))

lets_go <- 
vaf_ready %>%
  group_by(sample, gene) %>% 
  summarise(
    rna_vaf = sum(rna_vaf_est*rna_dp)/sum(rna_dp),
    rna_dp = sum(rna_dp),
    rna_dp_adj = sum(rna_dp * purple_af_mx),
    purple_af_mn = mean(purple_af_mx),
    vafs = n(),
    biallelic = paste0(unique(biallelic), collapse = ","), 
    afs = paste0(purple_af, collapse = ",")  
  ) %>% ungroup() 

fwrite(lets_go, "/home/jusset/vaf_clean/data/vaf_ready.txt")
