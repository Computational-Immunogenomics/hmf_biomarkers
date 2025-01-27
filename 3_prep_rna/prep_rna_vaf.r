source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

rna_vaf <- 
fread(paste0(TMP_DIR, "rna_vaf/rna_vaf.csv")) %>% 
 se(sampleId, chromosome, position, gene, rna_vaf, rna_dp)

SOM_DIR <- paste0(TMP_DIR, "somatic_exome_big/")

malig_exp <- data.frame()

system.time(
for( i in list.files(SOM_DIR)){
  print(i); flush.console()
  tmp <- 
  fread(paste0(SOM_DIR, i)) %>% 
  fi(purple_cn > .5, subclonal < .1, !is.na(purple_af)) %>% 
  se(sampleId, chromosome, position, purple_af, biallelic) %>% 
  ij(rna_vaf, by = c("sampleId" , "chromosome", "position")) %>% 
  mu(purple_af_mx = ifelse(biallelic, 1, purple_af), 
     rna_vaf_est = ifelse(rna_vaf/purple_af_mx > 1, 1, rna_vaf/purple_af_mx)) %>% 
  gb(sampleId, gene) %>% 
  su(
   malig_gene_pct = sum(rna_vaf_est*rna_dp)/sum(rna_dp),
  # rna_dp = sum(rna_dp), 
   rna_dp_adj = sum(rna_dp * purple_af_mx)
  # purple_af_mn = mean(purple_af_mx),
  # vafs = n(),
  # biallelic = paste0(unique(biallelic), collapse = ","), 
  # afs = paste0(purple_af, collapse = ",")
  ) %>% 
  ug()

  malig_exp <- rbind(malig_exp, tmp )
})

fwrite(malig_exp, paste0(TMP_DIR, "rna_vaf/malig_gene_exp.csv")) 
