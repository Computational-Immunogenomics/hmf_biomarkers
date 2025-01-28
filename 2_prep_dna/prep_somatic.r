source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

rna_vaf <- fread(paste0(TMP_DIR, "rna_vaf.csv"))

EXOME_DIR <- paste0(TMP_DIR, "somatic_exome_big/")

remove <- c("3_prime_UTR_variant","5_prime_UTR_variant","intron_variant","upstream_gene_variant","non_coding_transcript_exon_variant")

go <- data.frame()

for( i in list.files(EXOME_DIR)){
  fp <- paste0(EXOME_DIR, i)
  print(fp); flush.console()
  tmp <- 
  fread( fp ) %>% 
    fi(! annotation %in% remove) %>% 
    lj(rna_vaf, by = c("sampleId", "chromosome", "position", "REF", "ALT", "gene", "transcript"))
  go <- rbind(go, tmp)
}

fwrite(go, paste0(TMP_DIR, "somatic_exome.csv"))
