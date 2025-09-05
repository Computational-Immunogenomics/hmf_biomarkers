source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

linx_files <- c("breakend", "fusion", "svs", "vis_sv_data", "vis_fusion", "vis_protein_domain", "clusters", 
                "drivers", "links", "vis_copy_number", "vis_gene_exon", "vis_segments")

print("Collect SVs")

system.time(
for( j in linx_files){
  tmp <- list(); print(j); flush.console()
  for( i in list.files(SOM_DIR)){
    i_file <- get_linx_fp(i, type = j)
    if(file.exists(i_file)){ 
      tmp[[i]] <- reader(i_file, i) 
    }}
  fwrite( do.call("rbind", tmp), paste0(SV_DIR, j, ".csv"))    
}
)
