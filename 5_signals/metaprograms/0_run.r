source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

base <- fread(paste0(SHARE_DIR, "biomarkers_base.csv"))

cohorts <- fread("/mnt/petasan_immunocomp/datasets/hartwig/metadata/cohorts/cohorts_ready.csv")

ready <- 
base %>% 
 se( sampleId, 
     contains("cider_"), contains("clin_"), contains("cn_"), contains("driver"), contains("fusion_"), contains("gie_"), contains("lilac_"), 
     contains("purity_"), contains("sv_"), contains("teal_"), contains("viral_"), contains("bacterial_"), contains("metaprogram_")) %>% 
 drop_na(metaprogram_MP1_Cell_Cycle_G2M) %>% 
 lj(cohorts %>% se(sampleId, cohort), by = "sampleId")  %>% 
 se(where(~n_distinct(.) > 1)) %>% 
 se(where(~ !all(. %in% c(0, NA))))

go <- rbind(ready, ready %>% mu(cohort = "Pan-Cancer"))

metaprograms <- names(ready %>% se(contains("metaprogram_activity")))
features <- names(ready %>% se(-sampleId, -contains("metaprogram"), -cohort))
cohorts <- go %>% gb(cohort) %>% su(ct = n()) %>% fi(ct > 30) %>% ar(desc(ct)) %>% pu(cohort)

#go <- go %>% mu(across(any_of(features), scale))

results <- data.frame()
for(i in metaprograms){
 for( j in cohorts ) {
 run <- go %>% fi(cohort == j)
 if( j == "Pan-Cancer"){ 
  print(i); flush.console()   
  oo <- scanner(y = i, features, covariates = "+as.factor(cohort)", df = "run", mod = "lm") %>% mu(cohort = j)
 } else { 
  oo <- scanner(y = i, features, covariates = "", df = "run", mod = "lm") %>% mu(cohort = j)
 }
 results <- rbind(results, oo)
}}

fwrite(go , paste0(SHARE_DIR, "metaprogram_base_data.csv"))

fwrite(results %>% se(-lrt_pval, -data) , paste0(SHARE_DIR, "metaprogram_example_results.csv"))
