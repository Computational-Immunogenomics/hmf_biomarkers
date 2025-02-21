source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

annotation <- fread(paste0(dirname(getwd()), "/util/annotation/pathways.csv"))

annotation <- fread(paste0(dirname(getwd()), "/util/annotation/pathways.csv"))

drivers <- 
fread( paste0(TMP_DIR, "drivers.csv")) %>% 
  lj(annotation, by = "gene") %>% 
  mu(likelihood = ifelse(category == "TSG" & likelihoodMethod == "DISRUPTION", 1, driverLikelihood)) %>% 
  fi(likelihood > .8)  

drivers_ready <- 
drivers %>% 
  transmute(sampleId, gene = paste0("driver_", gene)) %>% 
  unique() %>% 
  mutate(driver = 1) %>% 
  spread(gene, driver)

drivers_ready[is.na(drivers_ready)] <- 0

pathways_ready <- 
drivers %>% 
 gb(sampleId, pathway) %>% 
 su(tot = n()) %>% 
 ug() %>% 
 tm(sampleId, pathway = paste0("drivers_pathway_", pathway), tot) %>% 
 unique() %>% 
 spread(pathway, tot)

pathways_ready[is.na(pathways_ready)] <- 0

total_drivers <- 
drivers %>% 
  gb(sampleId) %>% 
  su(drivers_total = n()) 

total_pathways <- 
drivers %>% 
 gb(sampleId) %>% 
 su(drivers_pathway_total = n_distinct(pathway))

together <- 
drivers_ready %>% 
 ij(pathways_ready, by = "sampleId") %>% 
 ij(total_drivers, by = "sampleId") %>% 
 ij(total_pathways, by = "sampleId")

fwrite(together, paste0(READY_DIR, "drivers_ready.csv"))
