source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "fisher.r"))

#library(ggh4x)
#library(patchwork)

ready <- readRDS(paste0(SHARE_DIR, "biomarkers_ready.Rds"))$ready 

cohorts <- 
fread("/mnt/bioinfnas2/immunocomp/shared_reference_data/cohorts/cohorts_ready.csv") %>% 
 se(sampleId, cohort) %>% 
 mu(cohort = ifelse( cohort %in% c("Colon", "Rectum"), "Colorectum", cohort))

go_treat <- 
rbind(ready %>% lj(cohorts, by = "sampleId"), ready %>% mu(cohort = "Pan-Cancer")) %>% 
 mu(cohortGo = paste0(cohort, " ## ", derived_treatmentName), group = "treatment")

go_mechanism <- 
rbind(ready %>% lj(cohorts, by = "sampleId"), ready %>% mu(cohort = "Pan-Cancer")) %>% 
 mu(cohortGo = paste0(cohort, " ## ", derived_treatmentMechanism), group = "mechanism")

go_cpi <- 
rbind(ready %>% lj(cohorts, by = "sampleId"), ready %>% mu(cohort = "Pan-Cancer")) %>% 
 fi(groupedTreatmentType %in% "Immunotherapy") %>% 
 mu(cohortGo = paste0(cohort, " ## ", "Immune Checkpoint Inhibitor"), group = "mechanism" )

go_type <- 
rbind(ready %>% lj(cohorts, by = "sampleId"), ready %>% mu(cohort = "Pan-Cancer")) %>% 
 mu(cohortGo = paste0(cohort, " ## ", groupedTreatmentType), group = "type" ) %>% 
 fi(groupedTreatmentType %in% c("Chemotherapy", "Immunotherapy", "Targeted therapy", "Hormonal therapy"))

all <- rbind(ready %>% lj(cohorts %>% se(sampleId, cohort), by = "sampleId")) %>% mu(cohortGo = "Pan-Cancer", group = "type" )

remove <- c("Pan-Cancer ## Anti-AR", "Pan-Cancer ## Folinic acid ## Platinum ## Pyrimidine (ant)agonist ## Topoisomerase inhibitor", 
           "Pan-Cancer ## Abiraterone", "Pan-Cancer ## Fluorouracil ## Irinotecan ## Leucovorin ## Oxaliplatin", 
           "Unknown primary (e.g. CUP) ## Chemotherapy", 
           "Pancreas PAAD ## Folinic acid ## Platinum ## Pyrimidine (ant)agonist ## Topoisomerase inhibitor", 
           "Prostate ## Hormonal therapy", 
           "Pan-Cancer ## Bevacizumab ## Capecitabine ## Oxaliplatin", 
           "Colorectum ## Anti-VEGF ## Platinum ## Pyrimidine (ant)agonist")

cohort_maps <- 
c("Pancreas PAAD ## Fluorouracil ## Irinotecan ## Leucovorin ## Oxaliplatin" = "Pancreas PAAD ## FOLFIRINOX",
  "Colorectum ## Bevacizumab ## Capecitabine ## Oxaliplatin" = "Colorectum ## CAPEOX + Bevacizumab")

go <- 
go_treat %>% 
 bind_rows(go_mechanism) %>% bind_rows(go_cpi) %>% bind_rows(go_type) %>% bind_rows(all) %>% 
 fi(!cohortGo %in% remove) %>% 
 mu(cohortGo = ifelse(cohortGo %in% names(cohort_maps), cohort_maps[cohortGo], cohortGo), 
    pan = grepl("Pan-Cancer", cohortGo))

fwrite(go, paste0(SHARE_DIR, "fisher_base.csv"))

min_patients <- 30; min_response <- 15

top_mechanisms <- 
go %>% 
 gb(cohortGo, group) %>% 
 su(ct = n(), no_dcb = sum(non_response), dcb = ct - no_dcb) %>% 
 fi(ct > min_patients, no_dcb >= min_response, dcb >= min_response) %>% 
 ug()

fwrite(top_mechanisms, paste0(SHARE_DIR, "top_mechanisms.csv"))
