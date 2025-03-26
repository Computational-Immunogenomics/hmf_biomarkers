source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

ready <- fread(paste0(SHARE_DIR, "biomarkers_ready.csv"))

cohorts <- fread("/mnt/bioinfnas2/immunocomp/shared_reference_data/cohorts/cohorts_ready.csv")

features <- 
names(ready %>% se( contains("cider_"), contains("clin_"), contains("cn_"), contains("driver_"), 
                    contains("fusion_"), contains("gie_"), contains("lilac_"), contains("neo_"), 
                    contains("purity_"), contains("rna_"), contains("sv_"), contains("teal_"), contains("viral_")))

go <- 
rbind(ready %>% lj(cohorts %>% se(sampleId, cohort), by = "sampleId"), 
      ready %>% mu(cohort = "Pan-Cancer")) %>% 
 mu(cohortGo = paste0(cohort, " ## ", derived_treatmentMechanism))

top_mechanisms <- 
go %>% 
 gb(cohortGo) %>% 
 su(ct = n(), no_dcb = sum(nrDcb), dcb = ct - no_dcb) %>% 
 fi(ct > 40, no_dcb > 15, dcb > 15) %>% 
 fi(cohortGo != "Pan-Cancer ## Anti-AR")

scanner <- function (y = "nrBor", features, covariates, df = "df", mod = "coxph") {
 out <- data.frame()
 for (f in features) {
  if (grepl("rna_", f) | grepl("cdr3_", f)) {
   tmp <- get_stats2(y = y, x = f, covariate = paste0("+purity_purity", covariates), data = df, model = mod)}
  else {tmp <- get_stats2(y = y, x = f, covariate = covariates, data = df, model = mod)} 
  if (is.data.frame(tmp)) out <- rbind(out, tmp)}
 out
}

cohorts_test <- top_mechanisms %>% pull(cohortGo)

mle_results <- data.frame()
for(i in cohorts_test){
 for( j in c("nrBor", "nrDcb")){
 print(i); flush.console()
 if(grepl("Pan-Cancer", i)){ covariates = "+cohort+ clin_age_cont" }
 else{ covariates = "+ clin_age_cont"}
 df <- go %>% fi(cohortGo == i)    
 tmp <- scanner(y = j, features, covariates = covariates, df = "df", mod = "bor")   
 mle_results <- rbind(mle_results, tmp %>% mu(cohort = i))
}}

fwrite(mle_results, paste0(SHARE_DIR, "0_run_lr.csv"))
