source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

base <- fread(paste0(SHARE_DIR, "biomarkers_base.csv"))

cohorts <- fread("/mnt/petasan_immunocomp/datasets/hartwig/metadata/cohorts/cohorts_ready.csv")

ready <- 
base %>% 
 se( sampleId, bestOverallResponse, osEvent, daysToOsEvent, 
     contains("cider_"), contains("clin_"), contains("cn_"), contains("driver"), contains("fusion_"), contains("gie_"), contains("lilac_"), 
     contains("purity"), contains("sv_"), contains("teal_"), contains("viral_"), contains("bacterial_"), contains("rna")) %>% 
 lj(cohorts %>% se(sampleId, cohort), by = "sampleId")  %>% 
 se(where(~n_distinct(.) > 1), where(~ !all(. %in% c(0, NA)))) %>% 
 drop_na(bestOverallResponse, daysToOsEvent, rna_geneset_KEGG_O_GLYCAN_BIOSYNTHESIS)

go <- rbind(ready, ready %>% mu(cohort = "Pan-Cancer"))

names(go) <- gsub("[^a-zA-Z0-9]", "_", names(go))

gps <- c("non_responders", "overall")
features <- names(go %>% se(-sampleId, -cohort, -bestOverallResponse, -osEvent, -daysToOsEvent))
cohorts <- go %>% gb(cohort) %>% su(ct = n()) %>% fi(ct > 15) %>% ar(desc(ct)) %>% pu(cohort)
covariates <- c("", "+ purity")

go <- go %>% mu(across(any_of(features), scale))

scanner <- function (y = "Surv(Y_os_days, Y_os_event)", features, covariates, 
    df = "df", mod = "coxph") {
    out <- data.frame()
    for (f in features) {
     tmp <- get_stats2(y = y, x = f, covariate = covariates, data = df, model = mod)
     if (is.data.frame(tmp)) out <- rbind(out, tmp)
    }
     out
}

results <- data.frame()
for(i in gps){
 for( k in covariates){
 
 if(i == "overall") { tmp <- go }
 else if ( i == "responders" ) { tmp <- go %>% fi(bestOverallResponse == 1) } 
 else { tmp <- go %>% fi(bestOverallResponse == 0) }
    
 for( j in cohorts ) {
 if(i == "overall") print(j); flush.console()  

 if(j == "Pan-Cancer") { run <- tmp %>% fi(cohort != j)}
 else {run <- tmp %>% fi(cohort == j)}

 if( i == "overall"){
   oo <- scanner(y = "Surv(daysToOsEvent, osEvent)", features, covariates=paste0("+ clin_age + bestOverallResponse", k), df = "run", mod = "coxph") 
 } else {
   oo <- scanner(y = "Surv(daysToOsEvent, osEvent)", features, covariates=paste0("+ clin_age", k), df = "run", mod = "coxph") 
 }
     
 results <- rbind(results , oo %>% mu(gp = i, cohort = j))
}}}

fwrite(go , paste0(SHARE_DIR, "prognostic_base_survival_data.csv"))

fwrite(results %>% se(-lrt_pval, -data) , paste0(SHARE_DIR, "prognostic_survival_results.csv"))
