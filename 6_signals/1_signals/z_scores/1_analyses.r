options(dplyr.summarise.inform = FALSE)
library(tidyverse)
library(data.table)
library(survival)
library(gridExtra)

source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(dirname(dirname(dirname(getwd()))),'/stats.r'))

system.time(go <- readRDS(paste0(UTIL_DIR, "zscores_tmp.csv")))
system.time(go2 <- readRDS(paste0(UTIL_DIR, "zscores2_tmp.csv")))

a <- ggplot( results, aes( x = est, y = -log10(pval_by), color = type)) + 
  geom_point() + 
  theme_classic() + 
  xlab("Log Hazard") + 
  ylab("-Log10 (BY Adjusted p-value)") + 
  ggtitle("Cox-ph: OS vs features (tissue, age, sex adjusted)") + 
  geom_hline(yintercept = -log10(.01)) + 
  facet_wrap(~type, ncol = 4) + 
  theme(legend.position = "none")

b <- ggplot( lms, aes( x = est, y = -log10(pval_by), color = type)) + 
  geom_point() + 
  theme_classic() + 
  xlab("Beta estimate") + 
  ylab("-Log10 (BY Adjusted p-value)") + 
  ggtitle("LM: Malignant_MP1 vs features (covariate adjusted)") + 
  geom_hline(yintercept = -log10(.01)) + 
  facet_wrap(~type, scales = "free", ncol = 4) + 
  theme(legend.position = "none")

#### 

options(repr.plot.width = 10, repr.plot.height = 4)
grid.arrange(a,b, ncol = 2)

df <- df %>% mutate( drivers_TP53_RB1 = (driver_TP53 > 0) + (driver_RB1 > 0))

a <- 
ggplot(df %>% drop_na(rna_mp_Malignant_MP1..Cell.Cycle...G2.M, drivers_TP53_RB1) %>% filter(clin_primaryTumorLocation2 != "Other"),
       aes( x =  as.factor(clin_primaryTumorLocation2), 
            y = rna_mp_Malignant_MP1..Cell.Cycle...G2.M,
            fill = as.factor(drivers_TP53_RB1))) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab("Malignant MP1 Cell Cycle") + 
  ggtitle("Malignant MP1 Cell Cycle vs Number of Drivers TP53 + RB1") + 
  theme(legend.position = "bottom")

b <- 
ggplot(df %>% drop_na(rna_mp_Malignant_MP1..Cell.Cycle...G2.M, drivers_TP53_RB1) %>% filter(clin_primaryTumorLocation2 != "Other"),
       aes( x =  as.factor(clin_primaryTumorLocation2), 
            y = rna_geneset_gene_set_prolif,
            fill = as.factor(drivers_TP53_RB1))) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab("RNA Proliferation Gene Set") + 
  ggtitle("RNA Proliferation Gene Set vs Number of Drivers TP53 + RB1") + 
  theme(legend.position = "bottom")

c <- 
ggplot(df %>% drop_na(rna_mp_Malignant_MP1..Cell.Cycle...G2.M, drivers_TP53_RB1) %>% filter(clin_primaryTumorLocation2 != "Other"),
       aes( x =  rna_geneset_gene_set_prolif, 
            y = rna_mp_Malignant_MP1..Cell.Cycle...G2.M,
            color = as.factor(clin_primaryTumorLocation2))) + 
  geom_point() + 
  theme_classic() + 
  ylab("Malignant MP1 Cell Cycle") + 
  xlab("RNA Proliferation Gene Set") + 
  ggtitle("Malignant MP1 Cell Cycle vs RNA Proliferation Gene Set") + 
  theme(legend.position = "bottom")

d <- 
ggplot(df %>% drop_na(rna_mp_Malignant_MP1..Cell.Cycle...G2.M, drivers_TP53_RB1) %>% filter(clin_primaryTumorLocation2 != "Other"),
       aes( x =  purity, 
            y = rna_mp_Malignant_MP1..Cell.Cycle...G2.M,
            color = as.factor(clin_primaryTumorLocation2))) + 
  geom_point() + 
  theme_classic() + 
  ylab("Malignant MP1 Cell Cycle") + 
  xlab("Scaled Purity") + 
  ggtitle("Malignant MP1 Cell Cycle vs Scaled Purity") + 
  theme(legend.position = "bottom")

options(repr.plot.width = 12, repr.plot.height = 10)
grid.arrange(a,b, c, d, ncol = 2)
