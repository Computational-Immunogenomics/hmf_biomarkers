source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

#rna_fusions <- fread(paste0(TMP_DIR, "isofox_fusions.csv"))

dna_fusions <- fread(paste0(TMP_DIR, "structural_variants/fusion.csv"))
#dna_fusions2 <- fread(paste0(TMP_DIR, "structural_variants/vis_fusion.csv"))

top_fusions <- dna_fusions %>% fi(geneStart != geneEnd) %>% gb(name) %>% su(ct = n()) %>% fi(ct > 50) %>% ar(desc(ct))

total_fusions <- 
dna_fusions %>% 
 gb(sampleId) %>% 
 su(total = n(), total_reported = sum(reported))

base <- 
dna_fusions %>% 
 fi(reported) %>% 
 tm(sampleId, name, fusion = 1 ) 

top_fusions <- 
base %>% 
 fi(name %in% (base %>% gb(name) %>% su(ct = n()) %>% fi(ct > 20) %>% ar(desc(ct)) %>% pu(name))) %>% 
 sp(name, fusion)

fusions_ready <- 
total_fusions %>% 
 lj(top_fusions, by = "sampleId") %>% 
 rename_with(~paste0("fusion_", .), .cols = -sampleId) %>% 
 mutate(across(everything(),~replace_na(., 0)))

fwrite(fusions_ready, paste0( READY_DIR, "fusions_dna_ready.csv"))
