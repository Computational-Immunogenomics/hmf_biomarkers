HELPER_DIR <- paste0(getwd(),'/helpers/')

source(paste0(HELPER_DIR,'map.r'))
source(paste0(HELPER_DIR,'shortcuts.r'))
source(paste0(HELPER_DIR,'helpers.r'))

base <- fread(paste0(SHARE_DIR, "biomarkers_base.csv"))

ready <- 
base %>% 
 se( sampleId, 
     cohort, 
     biopsyStructure, 
     contains("driver"), 
     contains("purity"), 
     contains("teal_")) %>% 
 se(where(~n_distinct(.) > 1)) %>% 
 se(where(~ !all(. %in% c(0, NA)))) %>% 
 mu(biopsy = ifelse(biopsyStructure %in% c("Liver", "Lymph node", "Bone", "Lung"), biopsyStructure, "Other")) %>%  
 mu(across(where(is.numeric), ~ replace_na(., median(., na.rm = TRUE)))) %>% 
 mu(sv_group = cut(purity_svTMB, breaks = 3, labels = c("Low", "Medium", "High")))

ready <- 
rbind(ready,
      ready %>% fi(sv_group == "Low") %>% mu(cohort = "Pan-Cancer: SV Burden Low"),
      ready %>% fi(sv_group == "Medium") %>% mu(cohort = "Pan-Cancer: SV Burden Med"),
      ready %>% fi(sv_group == "High") %>% mu(cohort = "Pan-Cancer: SV Burden High"))

non_epithelial <- 
ready %>%
 fi(grepl("Soft tissue", cohort) | cohort == "Glioblastoma" | grepl("NET", cohort) | grepl("Melanoma", cohort)) %>%
 mu(cohort = "Non-Epithelial")

epithelial <- 
ready %>%
 fi(!(grepl("Soft tissue", cohort) | cohort == "Glioblastoma" | grepl("NET", cohort) | grepl("Melanoma", cohort))) %>%
 mu(cohort = "Epithelial")

ready <- rbind(ready, non_epithelial, epithelial)

telomeres <- names(ready %>% se(contains("teal")))
features <- names(ready %>% se(-sampleId, -contains("teal"), -cohort, -biopsyStructure, -biopsy, -cohort, -sv_group))
cohorts <- c(ready %>% gb(cohort) %>% su(ct = n()) %>% fi(ct > 30) %>% ar(desc(ct)) %>% pu(cohort), "Pan-Cancer")
covariates <- c("", 
                " + purity_ploidy",
                "+ as.factor(biopsy)", 
                "+ as.factor(biopsy) + purity_ploidy", 
                "+ as.factor(biopsy) + purity", 
                "+ as.factor(biopsy) + purity + purity_ploidy",
                "+ as.factor(cohort) + as.factor(biopsy) + purity",
                "+ as.factor(cohort) + as.factor(biopsy) + purity + purity_ploidy"
               )

go <- ready %>% mu(across(any_of(features), scale))

results <- data.frame()
#for(i in c(telomeres, "purity_svTMB")){
for(i in c("teal_final_ratio", "purity_svTMB")){
 print(i); flush.console()     
 for( j in cohorts ) {
  print(j); flush.console()
  if(j == "Pan-Cancer"){ run <- go }
  else { run <- go %>% fi(cohort == j)}
  for( k in covariates){
    print(k); flush.console()   
    oo <- scanner(y = i, features, covariates = k, df = "run", mod = "lm")
 results <- rbind(results, oo %>% mu(cohort = j))
}}}

fwrite( results, paste0("data/0_run.csv"))
