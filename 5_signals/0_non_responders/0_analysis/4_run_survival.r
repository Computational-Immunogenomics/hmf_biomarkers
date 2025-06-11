source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

ready <- fread(paste0(SHARE_DIR, "biomarkers_ready.csv")) %>% fi(!is.na(durableClinicalBenefit))

cohorts <- fread("/mnt/bioinfnas2/immunocomp/shared_reference_data/cohorts/cohorts_ready.csv")

features <- 
names(ready %>% se( contains("cider_"), contains("clin_"), contains("cn_"), contains("driver_"), 
                    contains("fusion_"), contains("gie_"), contains("lilac_"), contains("neo_"), 
                    contains("purity_"), contains("rna_"), contains("sv_"), contains("teal_"), contains("viral_")))

base <- ready %>% lj(cohorts %>% se(sampleId, cohort), by = "sampleId")
go <- 
rbind(base, base %>% mu(cohort = paste0("Pan-Cancer ## ", cohort))) %>% 
 mu(cohortGo = paste0(cohort, " ## ", derived_treatmentMechanism))

top_mechanisms <- 
go %>% 
 gb(cohortGo) %>% 
 su(ct = n(), no_dcb = sum(nrDcb), dcb = ct - no_dcb) %>% 
 fi(ct > 40, no_dcb > 15, dcb > 15) %>% 
 fi(cohortGo != "Pan-Cancer ## Anti-AR")

cohorts_test <- top_mechanisms %>% pull(cohortGo)

mle_results <- data.frame()
for(i in cohorts_test){
 for( j in c("nrBor", "nrDcb")){
 print(i); flush.console()
 if(grepl("Pan-Cancer", i)){ covariates = "+ clin_age_cont" }
 else{ covariates = "+ clin_age_cont"}
 df <- go %>% fi(cohortGo == i)    
 tmp <- scanner(y = j, features, covariates = covariates, df = "df", mod = "bor")   
 mle_results <- rbind(mle_results, tmp %>% mu(cohort = i))
}}

fwrite(mle_results, paste0(SHARE_DIR, "1_run_lr.csv"))
