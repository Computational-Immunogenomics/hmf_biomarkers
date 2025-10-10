HELPER_DIR <- paste0(getwd(),'/helpers/')

source(paste0(HELPER_DIR,'map.r'))
source(paste0(HELPER_DIR,'shortcuts.r'))
source(paste0(HELPER_DIR,'helpers.r'))

base <- fread(paste0(SHARE_DIR, "biomarkers_base.csv"))

cohorts <- 
fread("/mnt/bioinfnas2/immunocomp/shared_reference_data/cohorts/cohorts_ready.csv") %>% 
 se(sampleId, cohort) %>% 
 mu(cohort = ifelse( cohort %in% c("Colon", "Rectum"), "Colorectum", cohort))

ready <- 
base %>% 
 se( sampleId, 
     biopsyStructure, 
     contains("driver"), 
     contains("fusion_"),
     contains("purity"), 
     contains("teal_"), 
     contains("viral_"), 
     contains("hotspot")) %>% 
 lj(cohorts %>% se(sampleId, cohort), by = "sampleId")  %>% 
 se(where(~n_distinct(.) > 1)) %>% 
 se(where(~ !all(. %in% c(0, NA)))) %>% 
 mu(across(where(is.numeric), ~ replace_na(., median(., na.rm = TRUE))),
    biopsy = ifelse(biopsyStructure %in% c("Liver", "Lymph node", "Bone", "Lung"), biopsyStructure, "Other")) 

telomeres <- names(ready %>% se(contains("teal")))
features <- names(ready %>% se(-sampleId, -contains("teal"), -cohort, -biopsyStructure, -biopsy))
cohorts <- c(ready %>% gb(cohort) %>% su(ct = n()) %>% fi(ct > 30) %>% ar(desc(ct)) %>% pu(cohort), "Pan-Cancer")
covariates <- c("", "+ as.factor(biopsy)", "+ as.factor(biopsy) + purity", "+ as.factor(cohort) + as.factor(biopsy) + purity")

go <- ready %>% mu(across(any_of(features), scale))

results <- data.frame()
for(i in telomeres){
 for( j in cohorts ) {
  if(j == "Pan-Cancer"){ run <- go }
  else { run <- go %>% fi(cohort == j)}
  for( k in covariates){
    print(i); flush.console()   
    #print(j); flush.console()
    oo <- scanner(y = i, features, covariates = k, df = "run", mod = "lm")
 results <- rbind(results, oo %>% mu(cohort = j))
}}}

fwrite( results, paste0("data/0_run.csv"))
