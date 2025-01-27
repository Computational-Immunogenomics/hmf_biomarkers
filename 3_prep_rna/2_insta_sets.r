options(repr.matrix.max.cols=50, repr.matrix.max.rows=100)
library(tidyverse)
library(data.table)
library(InstaPrism)

source(paste0(dirname(getwd()),'/map.r'))
gene_sets <- readRDS(paste0(REF_DIR, 'gene_sets_full.Rds'))

REF_DIR

names(gene_sets)

setwd(INST_DIR)

hires <- list()
for( i in list.files()){
  if( grepl("_ct_", i) & !grepl("_Qian_",i)){
     load(i)
     hires[[i]][['Malignant']] <- 
     t(reconstruct_Z_ct_initial(InstaPrism_obj = results_ct, cell.type.of.interest = "Malignant"))
}}

computer <- function( i, df ) {
  tmp <- data.frame( apply(log(data.frame(df) %>% select(any_of(gene_sets[[i]]))+1),1,mean) )
  colnames(tmp) <- i
  tmp %>% rownames_to_column("sampleId")
}

computed_sets <- list()
system.time(
for( j in names(hires)){
  tmp <- list()  
  for( i in names(gene_sets)){ 
    tmp[[i]] <- computer(i, hires[[j]][['Malignant']])
  }
  computed_sets[[j]] <- tmp %>% reduce(inner_join, by = "sampleId")
})

malignant_sets_ready <- do.call("bind_rows", computed_sets)
colnames(malignant_sets_ready) <- c("sampleId", paste0("rna_malig_geneset_", colnames(malignant_sets_ready)[-1]))

fwrite( malignant_sets_ready, paste0(READY_DIR, "isofox_malignant_genesets_ready.csv"))
