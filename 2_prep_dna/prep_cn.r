source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

cn <- fread(paste0(TMP_DIR, "cn_region.csv"))

cn_ready <- 
cn %>% 
  mu(cn_region = paste0("cn_chr", chromosome, "_", chromosomeBand)) %>% 
  se(-chromosome, -chromosomeBand) %>% 
  sp( cn_region, cn ) 

fwrite(cn_ready, paste0(READY_DIR, "cn_ready.csv"))

cn_ready_simple <- 
cn %>% 
 mu(pq = ifelse(grepl("p", chromosomeBand), "p", "q"), cn_region_pq = paste0("cn_simple_chr", chromosome, "_", pq)) %>% 
 se(-chromosome, -chromosomeBand) %>% 
 gb( sampleId, cn_region_pq ) %>% 
 su( cn = mean(cn) ) %>% 
 sp( cn_region_pq, cn ) 

fwrite(cn_ready_simple, paste0(READY_DIR, "cn_simple_ready.csv"))

#cn_genes <- fread(paste0(TMP_DIR, "cn_gene.csv"))
