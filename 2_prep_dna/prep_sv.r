source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

svs <- fread(paste0(TMP_DIR, "structural_variants/vis_sv_data.csv")) %>% fi(!Type %in% c("SGL", "INF")) 

svsTmb <- 
svs %>% 
 gb(sampleId, ChrStart, ChrEnd) %>% 
 su(svTmb = n()) %>%
 ug() %>% 
 mu(pos = ifelse( ChrStart != ChrEnd, paste0("svTmb_chr", ChrStart, "_", "chr", ChrEnd), paste0("svTmb_chr", ChrStart))) %>% 
 se(-ChrStart, -ChrEnd) %>% 
 se(sampleId, svTmb, pos) %>% 
 sp(pos, svTmb) 

svsTmb[is.na(svsTmb)] <- 0

ecdna <- 
svs %>% 
 gb(sampleId) %>% 
 su( svTmb_ecdna = sum(InDoubleMinute), ecdna = as.numeric((svTmb_ecdna > 0))) %>%
 ug()

svs_ready <- 
svsTmb %>% 
 lj(ecdna, by = "sampleId") %>% 
 rename_with(~ paste0("sv_", .), .cols = -sampleId)

fwrite(svs_ready, paste0(READY_DIR, "svs_ready.csv"))
