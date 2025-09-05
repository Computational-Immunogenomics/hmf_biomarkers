### Specify input directories for the Hartwig data ###
I_DIR = "/mnt/immunocompnas1/datasets/hartwig/"
SOM_DIR = paste0(I_DIR, "somatics/")
ISOFOX_DIR = paste0(I_DIR, "isofox/data_isofox/")
LILAC_DIR = paste0(I_DIR, "lilac/lilac_out/")
TEAL_DIR = paste0(I_DIR, "teal/")
NEO_DIR = paste0(I_DIR, "neo/data_neo/")
VIRAL_DIR = paste0(I_DIR, "viral_integration/data/datasets/")
CIDER_DIR = paste0(I_DIR, "cider/")
META_DIR = paste0(I_DIR, 'metadata/')

### Output directories ###
TMP_DIR = paste0(I_DIR, "biomarkers/database/")
SV_DIR = paste0(TMP_DIR, "structural_variants/")
READY_DIR = paste0(I_DIR, "biomarkers/ready/")
SHARE_DIR = paste0(I_DIR, "biomarkers/share/")

### Utility directories ### 
REF_DIR= paste0(I_DIR, "biomarkers/ref/")
HELP_DIR="/home/josephusset@vhio.org/biomarkers/helpers/"
FIG_DIR="/home/josephusset@vhio.org/biomarkers/util/figs/"

### File path mapper for samples
get_fp <- function(i, type = "purity"){
  if( type == "purity") { ac( paste0(SOM_DIR, i, "/purple/", i, ".purple.purity.tsv")) }
  else if (type == "drivers") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.driver.catalog.tsv")) }
  else if (type == "cnv_gene") { ac(paste0(SOM_DIR, i, "/purple/", i, ".purple.cnv.gene.tsv")) }
  else if (type == "isofox"){ ac(paste0(ISOFOX_DIR, i, "/", i, ".isf.gene_data.csv")) }
  else if (type == "isofox_fusion"){ ac(paste0(ISOFOX_DIR, i, "/", i, ".isf.fusions.csv"))}
  else if (type == "lilac"){ ac(paste0(LILAC_DIR, i, ".lilac.tsv")) }
  else if (type == "lilac_qc"){ ac(paste0(LILAC_DIR, i, ".lilac.qc.tsv")) }
  else if (type == "teal"){ ac(paste0(TEAL_DIR, i, ".teal.tellength.tsv")) }
  else if (type == "neo"){ ac(paste0(NEO_DIR, i, "/", i, ".neo.neoepitope.tsv.gz")) }
  else if (type == "neo_pep"){ ac(paste0(NEO_DIR, i,"/", i, ".neo.peptide_scores.tsv.gz")) }
  else if (type == "viral"){ ac(paste0(VIRAL_DIR, i, "/virus_interpreter/", i, ".virus.annotated.tsv"))}
}

get_linx_fp <- function(i, type = "breakend"){
  if (type == "breakend") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.breakend.tsv")) }
  else if (type == "fusion") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.fusion.tsv")) }
  else if (type == "svs") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.svs.tsv")) }
  else if (type == "vis_sv_data") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_sv_data.tsv")) }
  else if (type == "vis_fusion") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_fusion.tsv" )) }
  else if (type == "vis_protein_domain") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_protein_domain.tsv")) }
  else if (type == "clusters") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.clusters.tsv")) }
  else if (type == "drivers") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.drivers.tsv")) }
  else if (type == "links") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.links.tsv")) }
  else if (type == "vis_copy_number") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_copy_number.tsv")) }
  else if (type == "vis_gene_exon") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_gene_exon.tsv")) }
  else if (type == "vis_segments") { ac(paste0(SOM_DIR, i, "/linx/", i, ".linx.vis_segments.tsv")) }
}

### File reader 
reader <- function( i_file = "file_path", sample = "ACTN01020001T"){
    fread( i_file ) %>% mutate(sampleId = sample) 
}
