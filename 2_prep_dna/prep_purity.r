source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

purity_ready <- 
fread( paste0(TMP_DIR, "purities.csv") ) %>%
  tm(
   sampleId, 
   purity, 
   purity_diploidProportion = diploidProportion, 
   purity_ploidy = ploidy, 
   purity_polyclonalProportion = polyclonalProportion, 
   purity_WGD = ifelse(wholeGenomeDuplication, 1,0) , 
   purity_msIndelsPerMb = log(msIndelsPerMb+1), 
   purity_msStatus = ifelse(msStatus == "MSI", 1, 0), 
   purity_tml = log(tml+1),	
   purity_tmlStatus = ifelse(tmlStatus == "HIGH", 1, 0),	
   purity_tmbPerMb = log(tmbPerMb+1), 
   purity_tmbStatus = ifelse(tmbStatus == "HIGH", 1, 0), 
   purity_svTMB = log(svTumorMutationalBurden+1))

fwrite(purity_ready, paste0(READY_DIR, "purity_ready.csv"))
