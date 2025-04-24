source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

base <- fread(paste0(SHARE_DIR, "biomarkers_base.csv"))

cohorts <- fread("/mnt/petasan_immunocomp/datasets/hartwig/metadata/cohorts/cohorts_ready.csv")

ready <- 
base %>% 
 se( sampleId, biopsyStructure, 
     contains("cider_"), contains("clin_"), contains("cn_"), contains("driver"), contains("fusion_"), contains("gie_"), contains("lilac_"), 
     contains("purity"), contains("sv_"), contains("teal_"), contains("viral_"), contains("bacterial_"), contains("metaprogram_"), 
     contains("chord_"), contains("hotspot"), contains("neo_"), contains("signature")) %>% 
 drop_na(metaprogram_MP1_Cell_Cycle_G2M) %>% 
 lj(cohorts %>% se(sampleId, cohort), by = "sampleId")  %>% 
 se(where(~n_distinct(.) > 1)) %>% 
 se(where(~ !all(. %in% c(0, NA)))) %>% 
 mu(across(where(is.numeric), ~ replace_na(., median(., na.rm = TRUE))),
    biopsy = ifelse(biopsyStructure %in% c("Liver", "Lymph node", "Bone", "Lung"), biopsyStructure, "Other")) 

go <- rbind(ready, ready %>% mu(cohort = "Pan-Cancer"))

metaprograms <- names(ready %>% se(contains("metaprogram_activity")))
features <- names(ready %>% se(-sampleId, -contains("metaprogram"), -cohort, -biopsyStructure, -biopsy))
cohorts <- go %>% gb(cohort) %>% su(ct = n()) %>% fi(ct > 30) %>% ar(desc(ct)) %>% pu(cohort)
covariates <- c("", "+ as.factor(biopsy)", "+ as.factor(biopsy) + purity")
#covariates <- c("+ as.factor(biopsy)")

go <- go %>% mu(across(any_of(features), scale))

results <- data.frame()
for(i in metaprograms){
 for( j in cohorts ) {
  for( k in covariates){
   if( j == "Pan-Cancer"){ 
    run <- go %>% fi(cohort != j)   
    print(i); flush.console()   
    oo <- scanner(y = i, features, covariates = paste0("+as.factor(cohort)", k), df = "run", mod = "lm")
   } else {
    run <- go %>% fi(cohort == j)   
    oo <- scanner(y = i, features, covariates = k, df = "run", mod = "lm")
  }
 results <- rbind(results, oo %>% mu(cohort = j))
}}}

fwrite(go , paste0(SHARE_DIR, "metaprogram_base_data.csv"))

fwrite(results %>% se(-lrt_pval, -data, -model) , paste0(SHARE_DIR, "metaprogram_lm_results.csv"))
