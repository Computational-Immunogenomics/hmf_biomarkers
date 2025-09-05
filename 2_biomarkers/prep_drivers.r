source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

annotation <- fread(paste0(dirname(getwd()), "/util/annotation/pathways.csv"))

mapping <- list(
  "angiogenesis" = "ANGIOGENESIS",
  "apoptosis" = "APOPTOSIS",
  "cell cycle" = "CELL_CYCLE",
  "chromatin remodeling/epigenetic regulation" = "CHROMATIN",
  "dna damage response" = "DDR",
  "metabolism" = "METABOLISM",
  "hedgehog" = "HEDGEHOG",
  "hippo" = "HIPPO",
  "pd-1/pd-l1/ctla4" = "IMMUNE_EVASION",
  "jak/stat" = "JAK",
  "mapk/erk" = "MAPK",
  "myc" = "MYC",
  "nf-kb" = "NFKB",
  "notch" = "NOTCH",
  "nrf2" = "NRF2",
  "unknown" = "OTHER",
  "pi3k" = "PI3K",
  "dna damage repair" = "DDR",
  "tgfb" = "TGFB",
  "tp53" = "TP53",
  "wnt" = "WNT", 
  "rtk/ras" = "RTK_RAS", 
  "hormone receptor" = "HR", 
  "emt/metastasis" = "EMT_MET",
  "ubiquitination/proteostasis" = "UBQ_PROTEASIS",
  "spliceosome/splicing" = "SPLICING",
  "autophagy" = "AUTOPHAGY",
  "nf-κb" = "NFKB",
  'antigen presentation' = "ANTIGEN_PRESENTATION", 
  'tcr signaling pathway' = "TCR_SIGNALLING", 
  'interferon signaling' = "INTERFERON"    
)

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

total_drivers <- drivers %>% gb(sampleId) %>% su(drivers_total = n()) 
total_pathways <- drivers %>% gb(sampleId) %>% su(drivers_pathway_total = n_distinct(pathway))
total_immune_evasion <- drivers %>% gb(sampleId) %>% su(drivers_immune_evasion = sum(pathway == "IMMUNE_EVASION"))

together <- 
drivers_ready %>% 
 ij(pathways_ready, by = "sampleId") %>% 
 ij(total_drivers, by = "sampleId") %>% 
 ij(total_pathways, by = "sampleId") %>% 
 ij(total_immune_evasion, by = "sampleId")

fwrite(together, paste0(READY_DIR, "drivers_ready.csv"))
