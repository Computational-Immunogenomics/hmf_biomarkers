source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "fisher.r"))

ready <- readRDS(paste0(SHARE_DIR, "biomarkers_ready.rds"))$ready 

cohorts <- 
fread("/mnt/bioinfnas2/immunocomp/shared_reference_data/cohorts/cohorts_ready.csv") %>% 
 se(sampleId, cohort) %>% 
 mu(cohort = ifelse( cohort %in% c("Colon", "Rectum"), "Colorectum", cohort))

go_treat <- 
rbind(ready %>% lj(cohorts, by = "sampleId"), ready %>% mu(cohort = "Pan-Cancer")) %>% 
 mu(cohortGo = paste0(cohort, " / ", derived_treatmentName), group = "treatment")

go_mechanism <- 
rbind(ready %>% lj(cohorts, by = "sampleId"), ready %>% mu(cohort = "Pan-Cancer")) %>% 
 mu(cohortGo = paste0(cohort, " / ", derived_treatmentMechanism), group = "mechanism")

go_type <- 
rbind(ready %>% lj(cohorts, by = "sampleId"), ready %>% mu(cohort = "Pan-Cancer")) %>% 
 mu(cohortGo = paste0(cohort, " / ", groupedTreatmentType), group = "type" ) %>% 
 fi(groupedTreatmentType %in% c("Chemotherapy", "Immunotherapy", "Targeted therapy", "Hormonal therapy"))

all <- 
ready %>% 
 lj(cohorts %>% se(sampleId, cohort), by = "sampleId") %>% 
 mu(cohortGo = "Pan-Cancer", group = "type" )

remove <- 
c("Pan-Cancer / Anti-AR", 
  "Pan-Cancer / Folinic acid / Platinum / Pyrimidine (ant)agonist / Topoisomerase inhibitor", 
  "Pan-Cancer / Abiraterone", 
  "Pan-Cancer / Fluorouracil / Irinotecan / Leucovorin / Oxaliplatin", 
  "Unknown primary (e.g. CUP) / Chemotherapy", 
  "Pancreas PAAD / Folinic acid / Platinum / Pyrimidine (ant)agonist / Topoisomerase inhibitor", 
  "Prostate / Hormonal therapy", 
  "Pan-Cancer / Bevacizumab / Capecitabine / Oxaliplatin", 
  "Colorectum / Anti-VEGF / Platinum / Pyrimidine (ant)agonist", 
  "Kidney / Targeted therapy",
  "Bladder Urothelial / Immunotherapy",
  "Pan-Cancer / Anti-CTLA-4 / Anti-PD-1",
  "Bladder Urothelial / Anti-PD-1",
  "Lung NSCLC / Immunotherapy",
  "Pan-Cancer / Anti-CTLA-4 / Anti-PD-1",
  "Skin Melanoma / Anti-CTLA-4 / Anti-PD-1",
  "Soft tissue GST / Targeted therapy",
  "Breast ER+/HER- / Pyrimidine (ant)agonist",
  "Pan-Cancer")

cohort_maps <- 
c("Pancreas PAAD / Fluorouracil / Irinotecan / Leucovorin / Oxaliplatin" = "Pancreas PAAD / FOLFIRINOX",
  "Colorectum / Bevacizumab / Capecitabine / Oxaliplatin" = "Colorectum / CAPEOX + Bevacizumab")

go <- 
go_treat %>% 
 bind_rows(go_mechanism) %>% 
 bind_rows(go_type) %>% 
 bind_rows(all) %>% 
 mu(cohortGo = gsub("##", "/", cohortGo),
    cohortGo = ifelse(cohortGo %in% names(cohort_maps), cohort_maps[cohortGo], cohortGo), 
    pan = grepl("Pan-Cancer", cohortGo)) %>% 
 fi(!cohortGo %in% remove) 

min_patients <- 30; min_response <- 15

top_mechanisms <- 
go %>% 
 gb(cohortGo, group) %>% 
 su(ct = n(), no_dcb = sum(non_response), dcb = ct - no_dcb) %>% 
 fi(ct > min_patients, no_dcb >= min_response, dcb >= min_response) %>% 
 ug()

top_mechanisms %>% gb(group, grepl("Pan-Cancer", cohortGo)) %>% su(ct = n())

print("Summaries of cohort sizes overall")
top_mechanisms %>% su(mean(ct), mean(no_dcb), median(ct), median(no_dcb), min(ct), max(ct))

print("Summaries of cohort sizes by treatment group")
top_mechanisms %>% gb(group) %>% su(n(), mean(ct), mean(no_dcb), median(ct), median(no_dcb), min(ct), max(ct))

fwrite(top_mechanisms, paste0(SHARE_DIR, "top_mechanisms.csv"))
fwrite(go %>% fi(cohortGo %in% top_mechanisms$cohortGo), paste0(SHARE_DIR, "fisher_base.csv"))
fwrite(go %>% 
       fi(cohortGo %in% top_mechanisms$cohortGo) %>% 
       se(sampleId, primaryTumorLocation, cohortGo, durableClinicalBenefit), 
       paste0(SHARE_DIR, "vignettes_config.csv"))
