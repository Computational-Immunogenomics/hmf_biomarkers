source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

drivers <- fread( paste0(TMP_DIR, "drivers.csv"))

drivers_ready <- 
drivers %>% 
  mu(likelihood = ifelse(category == "TSG" & likelihoodMethod == "DISRUPTION", 1, likelihoodMethod)) %>% 
  fi(likelihood > 1) %>%
  transmute(sampleId, gene = paste0("driver_", gene)) %>% 
  unique() %>% 
  mutate(driver = 1) %>% 
  spread(gene, driver)

drivers_ready[is.na(drivers_ready)] <- 0

fwrite(drivers_ready, paste0(READY_DIR, "drivers_ready.csv"))
