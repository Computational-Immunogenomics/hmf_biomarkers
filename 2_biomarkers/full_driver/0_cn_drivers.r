source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

cn_gene <- 
fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/cn_gene.csv")

purities <- fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/purities.csv")

tmp <- 
cn_gene %>% 
 ij(purities %>% se(sampleId, ploidy), by = "sampleId") %>%
 fi(maxCopyNumber > (3*ploidy)) %>% 
 mu( method = ifelse(minCopyNumber > 3*ploidy, "AMP", "PARTIAL_AMP"))

fwrite(tmp, "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/drivers_full/amp.txt")

tmp <- 
cn_gene %>% 
 ij( purities %>% se(sampleId, ploidy), by = "sampleId" ) %>%
 fi( minCopyNumber < .5 ) %>% 
 mu( method = "DEL" )  

fwrite(tmp, "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/drivers_full/del.txt")
