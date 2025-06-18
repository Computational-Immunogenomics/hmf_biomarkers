source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "fisher.r"))

go <- fread(paste0(SHARE_DIR, "fisher_base.csv"))
cohorts <- fread(paste0(SHARE_DIR, "top_mechanisms.csv"))
categorical_features <- readRDS(paste0(SHARE_DIR, "biomarkers_ready.Rds"))$features

ra_ready <- 
go %>% 
 fi(cohortGo %in% (cohorts %>% pu(cohortGo))) %>% 
 se(cohortGo, non_response, any_of(categorical_features))

ra_go <- ra_formatter_and_test(ra_ready)

scanner <- function (y = "Surv(daysToPfsEvent, pfsEvent)", features, covariates, df = "df", mod = "coxph", cohort = "Pan-Cancer") {
    out <- data.frame()
    if(grepl("Pan-Cancer", cohort)) covariates = paste0(covariates, "+ as.factor(primaryTumorLocation)")
    for (f in features) {
        if (grepl("rna_", f)) {tmp <- get_stats2(y = y, x = f, covariate = paste0(covariates, "+ purity"), data = df, model = mod)} 
        else {tmp <- get_stats2(y = y, x = f, covariate = covariates, data = df, model = mod)}
        if (is.data.frame(tmp)) out <- rbind(out, tmp)
    }
    out
}

oo_survival <- data.frame()
for( i in cohorts$cohortGo ){
  print(i); flush.console();
  tmp <- go %>% fi(cohortGo == i)
  tmp_oo <- scanner(  feature = categorical_features, covariates = "", df = "tmp", cohort = i)
  if( nrow(tmp_oo) > 0 ){
   oo_survival <- rbind(oo_survival, tmp_oo %>% mu(cohortGo = i) %>% se(-lrt_pval, -data, -model))
  }
}

lets_go <- 
ra_go %>%
 lj(oo_survival %>% tm( cohortGo, feature = x, surv_est = est, surv_se = se, surv_pval = pval), 
    by = c("feature", "cohortGo")) %>% 
 lj(cohort, by = "cohortGo")

fwrite(lets_go, paste0(SHARE_DIR, "2_run_marginal_output.csv"))
