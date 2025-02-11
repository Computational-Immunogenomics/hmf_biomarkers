source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

library(dndscv)

cohorts <- fread(paste0(REF_DIR, "cohorts_ready.csv")) %>% se(sampleId, cohort)
top_cohorts <- c(cohorts %>% gb(cohort) %>% su(ct = n()) %>% fi(ct > 50) %>% pu(cohort), "Pan-Cancer")

mutations <- cohorts %>% lj(fread(paste0(TMP_DIR, "somatic_exome.csv")), by = "sampleId")

overall <- data.frame()
sel_cvs <- data.frame()

for( i in top_cohorts ) { 
  print(i); flush.console()
  if(i == "Pan-Cancer"){ 
    go <- mutations %>% tm(sampleID = sampleId, chr = chromosome, pos = position, ref = REF, mut = ALT)
  } else {
    go <- mutations %>% fi(cohort == i ) %>% tm(sampleID = sampleId, chr = chromosome, pos = position, ref = REF, mut = ALT)
  }
  dndsout = dndscv(go)  
  overall <- overall %>% bind_rows(dndsout$globaldnds %>% mu(cohort = i))
  sel_cvs <- sel_cvs %>% bind_rows(dndsout$sel_cv %>% mu(cohort = i))  
}

fwrite(overall, paste0(TMP_DIR, "dnds/overall.csv"))
fwrite(sel_cvs, paste0(TMP_DIR, "dnds/sel_cvs.csv"))
