source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

cider_dna <- fread(paste0(TMP_DIR, "cider_dna.csv")) 

cider_dna_ready <-
cider_dna %>%
  fi(filter == "PASS") %>% 
  se(sampleId, locus) %>% 
  gb(sampleId, locus) %>% 
  su(ct = n()) %>% 
  sp(locus, ct) %>% 
  mu(TOT = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% 
  ug()

cider_dna_ready[is.na(cider_dna_ready)] <- 0
colnames(cider_dna_ready) <- c("sampleId", paste0("cider_dna_", colnames(cider_dna_ready)[-1]))

fwrite(cider_dna_ready, paste0(READY_DIR, "cider_dna_ready.csv"))
