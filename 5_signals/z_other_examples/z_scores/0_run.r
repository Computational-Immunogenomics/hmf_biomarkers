options(dplyr.summarise.inform = FALSE)
library(tidyverse)
library(data.table)
library(survival)
library(gridExtra)

source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(dirname(dirname(dirname(getwd()))),'/stats.r'))

system.time(go <- readRDS(paste0(SHARE_DIR, "ready_ex.Rds")))

df <- go$data_ready
features <- go$features

results <- data.frame()
system.time(
for( i in features){
    results <- 
      rbind(results, 
            get_stats2( y = "Surv(Y_os_days, Y_os_event)", 
                        x = i, 
                        covariate = " + clin_primaryTumorLocation2 + clin_age + clin_sex", 
                        data = "df", 
                        model = "coxph"))
})

fwrite(results %>% arrange(pval), paste0(UTIL_DIR, "zscores_tmp.csv"))

results <- fread(paste0(UTIL_DIR, "zscores_tmp.csv"))

top_zscores <- results %>% filter(grepl("zscore", x)) %>% arrange(pval) %>% head(20) %>% pull(x)

lms <- data.frame()
for(j in top_zscores){
  print(j); flush.console()
  for( i in features){
    lms <- rbind(lms, 
                 get_stats2( y = j, 
                 x = i, 
                 covariate = " + clin_primaryTumorLocation2 + clin_age + clin_sex", 
                 data = "df", 
                 model = "lm"))
}}

lms <- lms %>% mutate(pval_by = p.adjust(pval, method = "BY")) 

fwrite(lms, paste0(UTIL_DIR, "zscores_tmp2.csv"))
