source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

cohorts <- fread("deconvolved_ids_Mets.csv")

malig_vaf <- 
fread(paste0(TMP_DIR, "rna_vaf/malig_gene_exp.csv")) %>% 
 fi(rna_dp_adj > 8)

iso <- 
fread(paste0(TMP_DIR, "isofox_adj_tmp.csv")) %>% 
 se(-GeneId) %>% 
 rename(gene = GeneName) %>% 
 ga(sampleId, tpm, -gene)

together <- 
base %>% 
 ij(go, by = c("sampleId", "gene")) %>% 
 mu( tpm, malig_tpm = tpm * malig_gene_pct, me_tpm = tpm - malig_tpm)

fwrite(together, paste0(TMP_DIR, "rna_vaf/malig_gene_exp_ready.csv"))

tmp <- 
cohorts %>% 
 lj(together, by = "sampleId") %>% 
 se(cancer_Type, sampleId, gene, tpm, malig_tpm, me_tpm) %>% 
 fi(tpm != 0) %>% 
 drop_na(tpm)

maker <- function( cohort = "Skin Melanoma") {
  me <- 
  tmp %>% 
   fi(cancer_Type == cohort) %>% 
   se(sampleId, gene, me_tpm) %>% 
   sp(sampleId, me_tpm) %>% 
   rename_with(~ paste0("non_malignant_TPM_", .), .cols = -gene)

  malig <- 
  tmp %>% 
   fi(cancer_Type == cohort) %>% 
   se(sampleId, gene, malig_tpm) %>% 
   sp(sampleId, malig_tpm) %>% 
   rename_with(~ paste0("malignant_TPM_", .), .cols = -gene)  

  malig %>% full_join(me, by = "gene") %>% mutate(across(everything(), ~ replace_na(., 0)))
  
}

vaf_based_references <- list()

for( i in unique(cohorts$cancer_Type)){
  print(i); flush.console()
  vaf_based_references[[i]] <- maker(i)
}

saveRDS(vaf_based_references, paste0(TMP_DIR, "rna_vaf/rna_vaf_based_deconvolution_references.Rds"))
