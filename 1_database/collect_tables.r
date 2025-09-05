source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

print("Collect Purities")

patients <- list.files(SOM_DIR)

tmp <- list()
system.time(
for( i in patients){
    i_file <- get_fp(i, type = "purity")
    if(file.exists(i_file)){
        tmp[[i]] <- reader(i_file, i)
    }
})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "purities.csv"))

print("Collect Drivers")

tmp <- list()
system.time(
for( i in patients){
    i_file <- get_fp(i, type = "drivers")
    if(file.exists(i_file)){
        tmp[[i]] <- reader(i_file, i)
    }
})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "drivers.csv"))

print("Collect Copy Number Region")

tmp <- list()
system.time(
for( i in patients){
  i_file <- get_fp(i, type = "cnv_gene")
  if(file.exists(i_file)){
    tmp[[i]] <- 
      reader(i_file, i) %>% 
        group_by(sampleId, chromosome, chromosomeBand) %>% 
        summarise(cn = mean(maxCopyNumber), cn_minor_allele = mean(minMinorAlleleCopyNumber)) %>% 
        ungroup()
}})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "cn_region.csv"))

print("Collect Copy Number Gene")

tmp <- list(); j <- 0
system.time(
for( i in patients){
  j <- j+1
  print(j); flush.console()
  i_file <- get_fp(i, type = "cnv_gene")
  if(file.exists(i_file)){
    tmp[[i]] <- 
      reader(i_file, i) %>% se(sampleId, chromosome, chromosomeBand, gene, minCopyNumber, maxCopyNumber, minMinorAlleleCopyNumber)
}})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "cn_gene.csv"))

print("Collect Lilac")

lilac_patients <- unique(unlist(lapply(list.files(LILAC_DIR), function(i) strsplit(i, ".lilac")[[1]][1])))

tmp <- list(); 
system.time(
for( i in lilac_patients){
  i_file <- get_fp(i, type = "lilac")
  if(file.exists(i_file)){ 
    tmp[[i]] <- reader( i_file, sample = i) 
  }  
})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "lilac.csv"))

tmp <- list()
system.time(
for( i in lilac_patients){
  i_file_qc <- get_fp(i, type = "lilac_qc")    
  if(file.exists(i_file_qc)){ 
      tmp[[i]] <- reader( i_file_qc, sample = i) 
  }   
})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "lilac_qc.csv"))

print("Collect Teal")

teal_patients <- unique(unlist(lapply(list.files(TEAL_DIR), function(i) strsplit(i, ".teal")[[1]][1])))

tmp <- list()
system.time(
for( i in teal_patients){
  i_file = get_fp(i, type = "teal")
  if(file.exists(i_file)){   
      tmp[[i]] <- reader( i_file, sample = i)
  }
})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "teal.csv"))

print("Collect Viral")

viral_patients <- list.files(VIRAL_DIR)

tmp <- list(); 
system.time(
for( i in viral_patients ){
  i_file <- get_fp(i, type = "viral")
  if(file.exists(i_file)){    
      tmp[[i]] <- reader( i_file, sample = i)
  }
})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "viral_add.csv"))
fwrite( fread( paste0(I_DIR, 'viral/data_viral_integration_filtered.tsv.gz')), paste0(TMP_DIR, "viral.csv"))

print("Collect CIDER")

fwrite(fread(paste0(CIDER_DIR, "hmf_cider_dna.vdj.tsv")), paste0(TMP_DIR, "cider_dna.csv"))
fwrite(fread(paste0(CIDER_DIR, "hmf_cider_rna.vdj.tsv")), paste0(TMP_DIR, "cider_rna.csv"))

print("Collect Neoepitopes")

neo_patients <- list.files(NEO_DIR)

tmp <- list(); 
system.time(
for( i in neo_patients ){
  i_file <- get_fp(i, type = "neo")
  if(file.exists(i_file)){    
      tmp[[i]] <- reader( i_file = get_fp(i, type = "neo"), sample = i)
  }
})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "neo.csv"))

tmp <- list(); 
system.time(
for( i in neo_patients ){
  i_file <- get_fp(i, type = "neo_pep")
  if(file.exists(i_file)){    
      tmp[[i]] <- reader( i_file, sample = i ) %>% arrange(desc(Score)) %>% head(100)
  }
})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "neo_pep.csv"))

print("Collect Isofox")

isofox_patients <- list.files(ISOFOX_DIR)

isofox_reader <- function( i_file, sample){
    fread( i_file ) %>% select(GeneId, GeneName, AdjTPM)  %>% setNames(c("GeneId","GeneName", sample)) 
}

tmp <- list()
system.time(
for( i in isofox_patients){
    i_file <- get_fp( i, type = "isofox" )
    if(file.exists(i_file)){
        tmp[[i]] <- isofox_reader(i_file, i) 
    }
})

fwrite(tmp %>% reduce(inner_join, by = c("GeneId", "GeneName")), paste0(TMP_DIR, "isofox_adj_tmp.csv"))

tmp <- list()
system.time(
for( i in isofox_patients){
    i_file <- get_fp( i, type = "isofox_fusion" )
    if(file.exists(i_file)){
        tmp[[i]] <- reader( i_file, sample = i)
    }
})

fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "isofox_fusions.csv"))
