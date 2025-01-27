options(dplyr.summarise.inform = FALSE)
library(tidyverse)
library(data.table)
library( readxl )

source(paste0(dirname(getwd()),'/map.r'))

ISO_DIR <- paste0(I_DIR, 'isofox/data_isofox/')
patients <- list.files(ISO_DIR)

iso <- fread(paste0(TMP_DIR, "isofox_adj_tmp.csv"))

iso_base <- 
log(data.frame(
  t(iso |> 
    select(-GeneId) |>
    column_to_rownames("GeneName"))
) + 1)

gene_sets <- readRDS(paste0(REF_DIR, 'gene_sets.Rds'))

cell_types <- c("B_cells", "Endothelial", "Epithelial", "Fibroblasts", "Macrophages","CD4", "CD8", "Malignant")

mps <- list()
for(i in cell_types){
    tmp <- read_excel(paste0(REF_DIR, "/41586_2023_6130_MOESM14_ESM.xlsx"), sheet = i)
    names(tmp) <- paste0("mp_", i, "_", names(tmp))
    gene_sets <- c(gene_sets, as.list(tmp))
}
saveRDS(gene_sets, paste0(REF_DIR, 'gene_sets_full.Rds'))

computer <- function( i, df ) {
  tmp <- data.frame( apply(df %>% select(any_of(gene_sets[[i]])),1,mean) )
  colnames(tmp) <- i
  tmp %>% rownames_to_column("sampleId")
}

computed_sets <- list()
system.time(
for( i in names(gene_sets)){ 
  computed_sets[[i]] <- computer(i, iso_base)
})

gene_sets_base <- computed_sets %>% reduce(inner_join, by = "sampleId")

isofox_ready <- iso_base
colnames(isofox_ready) <- paste0("rna_", colnames(iso_base))
isofox_ready <- isofox_ready %>% rownames_to_column("sampleId")

fwrite( isofox_ready, paste0(READY_DIR, "isofox_genes_ready.csv"))

gene_sets_ready <- gene_sets_base
colnames(gene_sets_ready) <- c("sampleId", paste0("rna_geneset_", colnames(gene_sets_base)[-1]))

fwrite( gene_sets_ready, paste0(READY_DIR, "isofox_genesets_ready.csv"))
