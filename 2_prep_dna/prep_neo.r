source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

neo <- fread( paste0(TMP_DIR, "neo.csv"))

neo_ready <- 
neo %>% 
 fi( VariantCopyNumber > .5, SubclonalLikelihood < .2) %>% 
 gb( sampleId ) %>% 
 su(neo_ct = n())

neo_pep <- fread( paste0(TMP_DIR, "neo_pep.csv"))

top_peps <- 
neo_pep %>% 
  gb(Peptide) %>% 
  su(ct = n()) %>% 
  ar(desc(ct)) %>% 
  fi(ct > 30 ) %>% 
  pu(Peptide)

neo_pep_ready <- 
neo_pep %>% 
 fi(Peptide %in% top_peps) %>% 
 gb(sampleId, Peptide) %>%
 mu( rk = row_number(desc(Score))) %>%
 fi( rk == 1 ) %>%
 se(sampleId, Peptide, Score) %>%
 sp(Peptide, Score)

neo_pep_ready[is.na(neo_pep_ready)] <- 0

colnames(neo_pep_ready) <- c("sampleId", paste0("neo_pep_", colnames(neo_pep_ready)[-1]))

fwrite(neo_ready, paste0(READY_DIR, "neo_ready.csv"))
fwrite(neo_pep_ready, paste0(READY_DIR, "neo_pep_ready.csv"))
