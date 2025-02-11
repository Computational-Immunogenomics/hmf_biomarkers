source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

lilac_qc <- fread( paste0(TMP_DIR, "lilac_qc.csv"))
lilac <- fread( paste0(TMP_DIR, "lilac.csv")) %>% filter(sampleId %in% (lilac_qc %>% filter(Status == "PASS") %>% pull(sampleId)))

lilac_base <-
lilac %>% 
 rw() %>% 
 mu(supertype = gsub("\\*", "", strsplit(Allele, ":")[[1]][1]), 
    receptor = strsplit(Allele, "\\*")[[1]][1], 
    loss = TumorCopyNumber < .5) %>% 
 ug() %>% 
 se(sampleId, receptor, Allele, supertype, TumorCopyNumber, loss)

hla_supertypes <- 
lilac_base %>% 
 gb(sampleId) %>% 
 su(hla_supertype = paste0(supertype, collapse = "__")) %>% 
 ug()

hla_supertypes_cn <- 
lilac_base %>% 
 gb(sampleId, supertype) %>%
 su(cn = mean(TumorCopyNumber)) %>% 
 sp(supertype, cn) %>% 
 mutate(across(everything(), ~ replace_na(., 0))) %>% 
 rename_with(~ paste0("hla_cn_", .), .cols = -sampleId) %>% 
 ug()

loh <- 
lilac_base %>% 
 gb(sampleId, receptor) %>% 
 su(tmp = sum(loss)) %>% 
 sp(receptor, tmp) %>% 
 mu(all = ifelse(A + B + C > 0, 1, 0)) %>% 
 ug() %>% 
 rename_with(~ paste0("hla_loh_", .), .cols = -sampleId)

zygosity <- 
lilac_base %>% 
 fi(!loss) %>% 
 gb(sampleId) %>% 
 su(hla_zygosity = n_distinct(supertype)) %>% 
 ug()

juntos <- 
hla_supertypes %>% 
 ij(hla_supertypes_cn, by = "sampleId") %>% 
 ij(loh, by = "sampleId") %>% 
 ij(zygosity, by = "sampleId") %>% 
 rename_with(~ paste0("lilac_", .), .cols = -sampleId)

fwrite(juntos, paste0(READY_DIR, "lilac_ready.csv"))
