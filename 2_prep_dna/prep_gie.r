source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

gie_mechanisms <- list(
    "HLA" = c("HLA-A", "HLA-B", "HLA-C"),
    "APM" = c("B2M", "CALR", "TAP1", "TAP2", "TAPBP", "NLRC5", "CIITA", "RFX5"),
    "IFG" = c("JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1"),
    "PDL1" = c("CD274"),
    "CD58" = c("CD58"),
    "SETDB1" = c("SETDB1")
)

gie_genes <- unname(unlist(gie_mechanisms))

ploidy <- fread(paste0(TMP_DIR, "purities.csv")) %>% se(sampleId, ploidy)
gie_cn <- fread(paste0(TMP_DIR, "cn_gene.csv")) %>% fi(gene %in% gie_genes)
gie_somatic <- fread(paste0(TMP_DIR, "somatic_exome.csv")) %>% fi(gene %in% gie_genes)
gie_svs <- fread(paste0(TMP_DIR, "structural_variants/breakend.csv")) %>% fi(gene %in% gie_genes, disruptive)

gie_cn_ready <- 
  gie_cn %>% 
    left_join(ploidy, by = "sampleId", relationship = "many-to-many") %>% 
    mutate(minMajoAlleleCopyNumber = minCopyNumber - minMinorAlleleCopyNumber) %>% 
     transmute(
          sampleId,
          gene,
          homozygous_del = ifelse(minCopyNumber < .5, TRUE, FALSE),
          loh =  ifelse( minMinorAlleleCopyNumber < 0.3 & minMajoAlleleCopyNumber > 0.7, TRUE, FALSE),
          amp = ifelse( minCopyNumber > 3 * ploidy, TRUE, FALSE),
          minCopyNumber,
          minMinorAlleleCopyNumber,
          minMajoAlleleCopyNumber, 
          ploidy) %>% 
    drop_na(ploidy) 

gie_somatic_ready <-
gie_somatic %>% 
  fi( subclonal	 < .5 ) %>% 
  mu( lof = ( grepl("stop_gained", annotation) |   grepl("frameshift", annotation) |  grepl("splice", annotation)), 
      bi_non_syn = ( biallelic & grepl("missense", annotation))) %>% 
  tm( sampleId, gene, lof, bi_non_syn) %>% 
  gb(sampleId, gene) %>% 
  su( lof = (sum(lof) > 0), bi_non_syn = (sum(bi_non_syn) > 0)) %>% 
  ug()

gie <-
  gie_cn_ready %>% 
    left_join(gie_somatic_ready, by = c("sampleId", "gene")) %>% 
    left_join(gie_svs %>% transmute(sampleId, gene, hd = TRUE ), by = c("sampleId", "gene")) %>% 
    select(sampleId, gene, homozygous_del, loh, amp, lof, bi_non_syn, hd) %>% 
    replace_na(list("lof" = FALSE, "bi_non_syn" = FALSE, hd = FALSE))

mapper <- list()
for( i in names(gie_mechanisms)){
    genes <- gie_mechanisms[[i]]
    for( j in genes ){
        mapper[[j]] <- i
    }
}   

gie$pathway <- unlist(lapply(gie$gene, function(i) mapper[[i]]))

gie_ready <-
  gie %>% 
    group_by(sampleId, pathway) %>% 
    summarise(h_del = sum(homozygous_del), 
              loh = sum(loh), 
              amp = sum(amp),
              lof = sum(lof),
              bi_non_syn = sum(bi_non_syn),
              hd = sum(hd)) %>% 
    mutate( gie_hla = (pathway == "HLA" & ( loh > 0 | h_del > 0 | lof > 0 | bi_non_syn > 0)), 
            gie_hla_lof = ( pathway == "HLA" & ( bi_non_syn > 0 | lof > 0 | hd > 0)),
            gie_hla_loh = ( pathway == "HLA" & ( loh > 0 | h_del > 0 )),
            gie_apm = ( pathway == "APM" & ( h_del > 0 | lof > 0 | bi_non_syn > 0 | hd > 0)),
            gie_ifg = ( pathway == "IFG" & ( h_del > 0 | lof > 0 | bi_non_syn > 0 | hd > 0)),
            gie_cd58 = ( pathway == "CD58" & ( h_del > 0 | lof > 0 | bi_non_syn > 0 | hd > 0)),
            gie_pdl1 = ( pathway == "PDL1" & ( amp )), 
            gie_setdb1 = ( pathway == "SETDB1" & ( amp )) ) %>% 
   group_by(sampleId) %>% 
   summarise( gie_hla = sum(gie_hla), 
              gie_hla_lof = sum(gie_hla_lof),
              gie_hla_loh = sum(gie_hla_loh),
              gie_apm = sum(gie_apm),
              gie_ifg = sum(gie_ifg), 
              gie_pdl1 = sum(gie_pdl1), 
              gie_cd58 = sum(gie_cd58),
              gie_setdb1 = sum(gie_setdb1))

gie_ready$gie <- ifelse(apply( gie_ready %>% select(contains("gie")), 1, sum) > 0, 1, 0)

fwrite( gie_ready, file = paste0( TMP_DIR, "gie.csv") )
