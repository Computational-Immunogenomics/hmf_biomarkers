source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

library(dndscv)

mutations <- 
fread(paste0(TMP_DIR, "somatic_exome.csv")) %>% 
  tm(sampleID = sampleId, chr = chromosome, pos = position, ref = REF, mut = ALT)

dndsout = dndscv(mutations)

saveRDS(dndsout, paste0(TMP_DIR, "dnds/dnds.Rds"))
