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

fwrite(ra_go %>% lj(cohorts, by = "cohortGo"), paste0(SHARE_DIR, "1_run_fishers_exact.csv"))
