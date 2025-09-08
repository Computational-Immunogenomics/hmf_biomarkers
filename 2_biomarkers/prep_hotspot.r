source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

somatic <- fread(paste0(TMP_DIR, "somatic_exome.csv")) 

threshold <- (1/200) * length(unique(somatic %>% pu(sampleId)))

non_coding <- c("upstream_gene_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant")

hot <- 
somatic %>% 
 fi(annotation != "synonymous_variant") %>% 
 gb(gene, annotation, chromosome, REF, ALT, position) %>% 
 su(ct = n()) %>% fi( ct >= threshold/5 ) %>% 
 gb(annotation) %>% mu(rk = row_number(desc(ct))) %>% 
 fi((!annotation %in% non_coding) | (annotation %in% non_coding & rk <= 20)) %>% 
 ug()

names_map <- 
c("hotspot_KRAS_missense_variant_chr12_refC_altT_pos25398284" = "hotspot_KRAS_G12D",
  "hotspot_KRAS_missense_variant_chr12_refC_altA_pos25398284" = "hotspot_KRAS_G12V",
  "hotspot_KRAS_missense_variant_chr12_refC_altG_pos25398284" = "hotspot_KRAS_G12A",
  "hotspot_KRAS_missense_variant_chr12_refC_altA_pos25398285" = "hotspot_KRAS_G12C",
  "hotspot_KRAS_missense_variant_chr12_refC_altG_pos25398285" = "hotspot_KRAS_G12G",
  "hotspot_KRAS_missense_variant_chr12_refC_altT_pos25398281" = "hotspot_KRAS_G13D", 
  "hotspot_BRAF_missense_variant_chr7_refA_altT_pos140453136" = "hotspot_BRAF_V600E",
  "hotspot_PIK3CA_missense_variant_chr3_refA_altG_pos178952085" = "hotspot_PIK3CA_H1047R",
  "hotspot_PIK3CA_missense_variant_chr3_refG_altA_pos178936091" = "hotspot_PIK3CA_E545K",
  "hotspot_PIK3CA_missense_variant_chr3_refG_altA_pos178936082" = "hotspot_PIK3CA_E542K",
  "hotspot_PIK3CA_missense_variant_chr3_refG_altA_pos178936082" = "hotspot_PIK3CA_E542K",
  "hotspot_TERT_upstream_gene_variant_chr5_refG_altA_pos1295228" = "hotspot_TERT_C228T",
  "hotspot_TERT_upstream_gene_variant_chr5_refG_altA_pos1295250" = "hotspot_TERT_C250T")
mapper <- function(i) if( i %in% names(names_map)){ names_map[[i]] } else{i}

maker <- function( hotspot_df, column = "position") {
 somatic %>% 
  ij( hotspot_df %>% se(gene, annotation, REF, ALT, position), by = c("gene", "position", "annotation", "REF", "ALT")) %>% 
  se( sampleId, gene, annotation, chromosome, REF, ALT, position ) %>% 
  mu( hotspot = paste0( "hotspot_", gene, "_", annotation, "_chr", chromosome, "_ref", REF, "_alt", ALT, "_pos", position ), ct = 1) %>% 
  se(-gene, -position, -annotation, -chromosome, -REF, -ALT) %>% 
  unique() %>% 
  rw() %>% mu( hotspot = mapper(hotspot)) %>% ug() %>% 
  sp(hotspot, ct)  
}

hotspots <- maker(hot)

hotspots_ready <- 
hotspots %>% 
 mu(across(everything(), ~ replace_na(., 0))) %>%
 mutate(hotspot_KRAS_G12 = as.numeric(rowSums(across(contains("KRAS_G12"))) > 0))

fwrite(hotspots_ready, paste0(READY_DIR, "hotspots_ready.csv"))
