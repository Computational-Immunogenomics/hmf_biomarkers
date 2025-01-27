source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

SOM_DIR <- paste0(I_DIR, "somatics/")

get_linx_fp <- function(i, type = "breakend"){
	if (type == "breakend") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.breakend.tsv")) }
	else if (type == "fusion") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.fusion.tsv")) }
	else if (type == "svs") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.svs.tsv")) }
	else if (type == "vis_sv_data") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_sv_data.tsv")) }
	else if (type == "vis_fusion") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_fusion.tsv" )) }
	else if (type == "vis_protein_domain") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_protein_domain.tsv")) }
	else if (type == "clusters") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.clusters.tsv")) }
	else if (type == "drivers") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.drivers.tsv")) }
	else if (type == "links") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.links.tsv")) }
	else if (type == "vis_copy_number") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_copy_number.tsv")) }
	else if (type == "vis_gene_exon") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_gene_exon.tsv")) }
	else if (type == "vis_segments") { as.character(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_segments.tsv")) }
}

linx_files <- c("breakend", "fusion", "svs", "vis_sv_data", "vis_fusion", "vis_protein_domain", "clusters", 
                "drivers", "links", "vis_copy_number", "vis_gene_exon", "vis_segments")

print("Collect SVs")

SOM_DIR <- paste0(I_DIR, "somatics/")
patients <- list.files(SOM_DIR)

system.time(
for( j in linx_files){
  tmp <- list()  
  print(j)
  flush.console()
  for( i in patients){
    i_file <- get_linx_fp(i, type = j)
    if(file.exists(i_file)){ 
      tmp[[i]] <- reader(i_file, i) 
    }}
  fwrite( do.call("rbind", tmp), paste0(TMP_DIR, "structural_variants/", j, ".csv"))    
}
)
