source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

teal_ready <- 
fread("/mnt/petasan_immunocomp/datasets/hartwig/teal/teal_summary.csv") %>% 
 tm( sampleId, 
     teal_ref = germlineTelomereLength, 
     teal_tumor = somaticTelomereLength, 
     teal_ratio = log2(teal_tumor/teal_ref+1) ) %>%
 se(-teal_ref)

fwrite(teal_ready, paste0(READY_DIR, "teal_ready.csv"))
