source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
library( readxl )

patients <- list.files(ISOFOX_DIR)

iso <- fread(paste0(TMP_DIR, "isofox_adj_tmp.csv"))

iso_base <- log(data.frame(t(iso %>% select(-GeneId) %>% column_to_rownames("GeneName"))) + 1)

mrp_genes <- c("ABCC1", "ABCC2", "ABCC3", "ABCC4", "ABCC5", "ABCC6", "ABCC10", "ABCC11", "ABCC12")
more_markers <- c("ESR1", "ERBB2", "CD274", "MGMT", "BRCA1", "BRCA2", "EGFR", "ALK", "BCL2", "AR", "TOP2A", "TYMS", "ERCC1", "MET", "KRAS")
marker_genes <- c(mrp_genes, more_markers)

biomarker_genes <-
iso_base %>% 
 se(any_of(marker_genes)) %>% 
 rownames_to_column("sampleId")

colnames(biomarker_genes) <- c("sampleId", paste0("rna_marker_", colnames(biomarker_genes)[-1]))

fwrite( biomarker_genes, paste0(READY_DIR, "biomarker_genes_ready.csv"))

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

fwrite( isofox_ready, paste0(READY_DIR, "isofox_genes_ready.csv") )

gene_sets_ready <- gene_sets_base
colnames(gene_sets_ready) <- c("sampleId", paste0("rna_geneset_", colnames(gene_sets_base)[-1]))

fwrite( gene_sets_ready, paste0(READY_DIR, "isofox_genesets_ready.csv"))
