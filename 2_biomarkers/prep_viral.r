source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

viral <- fread( paste0(TMP_DIR, "viral.csv") )

viral_reported <-
viral %>% 
 fi(reported) %>% 
 se(sampleId, interpretation, integrations) %>% 
 sp(interpretation, integrations) 

viral_reported[is.na(viral_reported)] <- 0

colnames(viral_reported) <- c("sampleId", paste0("viral_", colnames(viral_reported)[-1]))

fwrite(viral_reported, paste0(READY_DIR, "viral_ready.csv"))
