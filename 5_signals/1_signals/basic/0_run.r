options(dplyr.summarise.inform = FALSE)
library(tidyverse)
library(data.table)
library(survival)
library(gridExtra)

source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(dirname(dirname(dirname(getwd()))),'/stats.r'))

system.time(go <- readRDS(paste0(SHARE_DIR, "ready_ex.Rds")))

all <- go$data_ready
features <- go$features

cohorts <- list()
top_cohorts <- unique(all$location)
cohorts[['pan']] <- all
for( i in unique(all$location) ) cohorts[[i]] <- all %>% filter(location == i)

os_out <- data.frame()
for( c in names(cohorts)){
  df <- cohorts[[c]]
  if( c == "pan"){ os <- scanner("Surv(Y_os_days, Y_os_event)", features, "+ as.factor(location) + clin_age + clin_sex", "df", "coxph")} 
  else { os <- scanner("Surv(Y_os_days, Y_os_event)", features, "+ clin_age + clin_sex", "df", "coxph")}
  os_out <- rbind(os_out, os %>% mutate(cohort = c ))
}

bor_out <- data.frame()
for( c in names(cohorts)){
  df <- cohorts[[c]]
  if( c == "pan"){ bor <- scanner("Y_bor", features, "+ as.factor(location) + clin_age + clin_sex", "df", "bor")} 
  else { bor <- scanner("Y_bor", features, "+ clin_age + clin_sex", "df", "bor")}
  bor_out <- rbind(bor_out, bor %>% mutate(cohort = c ))
}

out <- rbind(os_out, bor_out) %>% mutate(pval_by = p.adjust(pval, method = "BY")) 

fwrite(out, paste0(UTIL_DIR, "go.csv"))
