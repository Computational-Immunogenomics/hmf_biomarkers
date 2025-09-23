source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
library( readxl )

patients <- list.files(ISOFOX_DIR)

iso <- fread(paste0(TMP_DIR, "isofox_adj_tmp.csv"))

iso_base <- log(data.frame(t(iso %>% select(-GeneId) %>% column_to_rownames("GeneName"))) + 1)

gene_sets <- readRDS(paste0(REF_DIR, 'gene_sets.Rds'))
gene_sets <- gene_sets[-grep("gene_set|vhio|battle|rand|character", names(gene_sets))]

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
