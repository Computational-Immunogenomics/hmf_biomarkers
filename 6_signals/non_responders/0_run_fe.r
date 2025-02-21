source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

ready <- fread(paste0(SHARE_DIR, "biomarkers_ready.csv"))

cohorts <- fread("/mnt/bioinfnas2/immunocomp/shared_reference_data/cohorts/cohorts_ready.csv")

categorical_features <- 
names(ready %>% se( contains("cider_"), contains("clin_"), contains("cn_"), contains("driver_"), 
                    contains("fusion_"), contains("gie_"), contains("lilac_"), contains("neo_"), 
                    contains("purity_"), contains("rna_"), contains("sv_"), contains("teal_"), contains("viral_"), 
                    -contains("_cont")))

top_mechanisms <- 
ready %>% 
 gb(primaryTumorLocation, derived_treatmentMechanism) %>% 
 su(ct = n(), no_dcb = sum(nrDcb), dcb = ct - no_dcb) %>% 
 fi(ct > 40, no_dcb > 15, dcb > 15) %>% 
 mu(cohort = "Pan-Cancer")

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

base <- 
go %>% 
 fi(cohortGo %in% (top_mechanisms %>% pu(cohortGo))) %>% 
 se(cohortGo, non_response = nrDcb, any_of(categorical_features)) %>% 
 ga(feature, event, -cohortGo, -non_response) %>% 
 drop_na(event) %>% 
 gb(cohortGo, feature, non_response, event) %>% 
 su(tot = n(), .groups = "drop") %>% 
 pivot_wider(names_from = c(event, non_response),  values_from = tot)

base[is.na(base)] <- 0

base <- base %>% mu(events = `1_0` + `1_1`) %>% fi( events > 5 )

ra_fisher <- function(a,b,c,d){
 fisher.test(matrix(c(a,b,c, d), ncol = 2))$p.value
}

ra_go <- 
base %>% 
 rw() %>% 
 mu(fisher_pval = ra_fisher(`0_0`, `0_1`, `1_0`, `1_1`)) %>% 
 ug() %>% 
 se(cohortGo, feature, `0_0`, `0_1`, `1_0`, `1_1`, events, fisher_pval)

ra_go %>% fi(events  == `1_1`) %>% ar(desc(events)) %>% fi(grepl("hrd", feature))

ra_ready <- 
ra_go %>% 
 ar(fisher_pval) %>% 
 rename("r_ne" = `0_0`, "nr_ne" = `1_0`, "r_e" = `0_1`, "nr_e" = `1_1`) %>%
 mu( tot_e = r_e + nr_e, 
     tot_ne = r_ne + nr_ne,
     tot_nr = nr_e + nr_ne, 
     tot_r = r_e + r_ne, 
     tot = tot_nr + tot_r, 
     pr_nr_given_e = nr_e/tot_e,
     pr_nr_overall = tot_nr/tot) %>% 
 se(cohortGo, feature, fisher_pval, 
    nr_e, r_e, nr_ne, r_ne, 
    tot_e, tot_ne, tot_nr, tot_r, tot, 
    pr_nr_overall, pr_nr_given_e) 

#dim(ra_ready)

fwrite(ra_ready, paste0(SHARE_DIR, "0_run_fe.csv"))
