source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

chord <- fread(paste0(I_DIR, "chord/chord.csv"))

chord_ready <-
chord %>% 
 mu(hrStatus = ifelse(hrStatus == "CANNOT_BE_DETERMINED", NA, hrStatus), 
    hrStatus = ifelse(hrStatus == "HR_DEFICIENT", 1, 0)) %>% 
 rename(chord_BRCA1 = BRCA1, chord_BRCA2 = BRCA2, chord_hrd = hrd, chord_hrStatus = hrStatus) %>% 
 se(sampleId, contains("chord"))

fwrite(chord_ready, paste0(READY_DIR, "chord_ready.csv"))
