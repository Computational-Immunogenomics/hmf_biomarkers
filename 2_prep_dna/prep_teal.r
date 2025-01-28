source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

teal_ready <- 
fread( paste0(TMP_DIR, "teal.csv")) %>% 
  mu( sampleId = gsub("R", "T", sampleId)) %>% 
  se( sampleId, type, finalTelomereLength) %>% 
  sp( type, finalTelomereLength) %>% 
  tm( sampleId, 
      teal_ref = ref, 
      teal_tumor = ifelse(tumor <= 0, 0, tumor), 
      teal_ratio = log2(teal_tumor/teal_ref+1) )

fwrite(teal_ready, paste0(READY_DIR, "teal_ready.csv"))
