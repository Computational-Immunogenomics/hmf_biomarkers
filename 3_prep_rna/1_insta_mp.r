options(repr.matrix.max.cols=50, repr.matrix.max.rows=100)
library(tidyverse)
library(data.table)

getwd()

source(paste0(dirname(getwd()),'/map.r'))

wd <- "/mnt/bioinfnas2/immunocomp/manuel/tme/metaPrograms/deconvolution_output_merged/"
zscores <- "/mnt/bioinfnas2/immunocomp/manuel/tme/metaPrograms/3_zscores_by_cellType/zscore_df.csv"

mps <- list()
for( i in list.files()){
  if( grepl("_mp_", i) & !grepl("_Qian_",i)){
     load(i)
     mps[[i]] <- data.frame(t(results_mp@Post.ini.cs@theta))
}}

together <- do.call("bind_rows", mps)
names(together) <- paste0( "rna_mp_", names(together))
mp_ready <- together %>% rownames_to_column("sampleId")

mp_zscores <- 
fread(zscores) %>% 
  mutate( gp = paste0(ct, "_", meta_program)) %>% 
  filter(sampleId != "TCGA_SAMPLE_433") %>% 
  select(-study, -ct, -tissue, -meta_program) %>% 
  spread(gp, zscore)
colnames( mp_zscores ) <- c("sampleId", paste0("rna_mp_zscore_", colnames(mp_zscores )[-1]))

cts <- list()
for( i in list.files()){
  if( grepl("_ct_", i) & !grepl("_Qian_",i)){
     load(i)
     cts[[i]] <- data.frame(t(results_ct@Post.ini.cs@theta))
}}

together <- do.call("bind_rows", cts)
names(together) <- paste0( "rna_ct_", names(together))
ct_ready <- together %>% rownames_to_column("sampleId")

fwrite( mp_ready, paste0(READY_DIR, "isofox_insta_mp_ready.csv"))
fwrite( mp_zscores, paste0(READY_DIR, "isofox_insta_mp_zscores_ready.csv"))
fwrite( ct_ready, paste0(READY_DIR, "isofox_insta_prop_ready.csv"))
