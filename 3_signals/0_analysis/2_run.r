source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "fisher.r"))

go <- fread(paste0(SHARE_DIR, "fisher_base.csv"))
cohorts <- fread(paste0(SHARE_DIR, "top_mechanisms.csv"))
categorical_features <- readRDS(paste0(SHARE_DIR, "biomarkers_ready.Rds"))$features

treatment_mechanism_map <- 
fread(paste0(SHARE_DIR, "treatment_mechanism_map.csv")) %>% 
 mu(derived_treatmentName = gsub( "##","/", derived_treatmentName),
    derived_treatmentMechanism = gsub( "##","/", derived_treatmentMechanism))

ra_ready <- 
go %>% 
 fi(cohortGo %in% (cohorts %>% pu(cohortGo))) %>% 
 se(cohortGo, non_response, any_of(categorical_features))

print("Run fisher's exact tests")
system.time(
 ra_go <- ra_formatter_and_test(ra_ready)
)

print("Add PFS survival analyses")
oo_survival <- data.frame()
system.time(
for( i in cohorts$cohortGo ){
  print(i); flush.console();
  tmp <- go %>% fi(cohortGo == i)
  tmp_oo <- scanner_non_responders(  feature = categorical_features, covariates = "", df = "tmp", cohort = i)
  if( nrow(tmp_oo) > 0 ){
   oo_survival <- rbind(oo_survival, tmp_oo %>% mu(cohortGo = i) %>% se(-lrt_pval, -data, -model))
  }
})

lets_go <- 
ra_go %>% 
 lj(oo_survival %>% 
     tm( cohortGo, feature = x, surv_est = est, surv_se = se, surv_pval = pval), 
     by = c("feature", "cohortGo"))

fwrite(lets_go, paste0(SHARE_DIR, "2_run_marginal_output.csv"))
