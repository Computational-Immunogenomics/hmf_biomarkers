I_DIR = "/mnt/petasan_immunocomp/datasets/hartwig/"
O_DIR = "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/"
TMP_DIR = "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/tmp/"
READY_DIR= "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/ready/"
SHARE_DIR = "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/share/"
REF_DIR= "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/ref/"
INST_DIR="/mnt/bioinfnas2/immunocomp/manuel/tme/metaPrograms/deconvolution_output_merged/"
UTIL_DIR="/mnt/bioinfnas/immunocomp/jusset/biomarkers/util/tmp/" 

### File path mapper for samples
get_fp <- function(i, type = "purity"){
  if( type == "purity") { as.character(paste0(I_DIR, 'somatic/', i, "/purple/", i, ".purple.purity.tsv")) }
  else if (type == "drivers") { as.character(paste0(I_DIR, 'somatic/', i, "/linx/", i, ".linx.driver.catalog.tsv")) }
  else if (type == "cnv_gene") { as.character(paste0(I_DIR, 'somatic/', i, "/purple/", i, ".purple.cnv.gene.tsv")) }
  else if (type == "isofox"){ as.character(paste0(I_DIR, 'isofox/data_isofox/', i, "/", i, ".isf.gene_data.csv")) }
  else if (type == "lilac"){ as.character(paste0(I_DIR, 'lilac/lilac_out/', i, ".lilac.tsv")) }
  else if (type == "lilac_qc"){ as.character(paste0(I_DIR, 'lilac/lilac_out/', i, ".lilac.qc.tsv")) }
  else if (type == "teal"){ as.character(paste0(I_DIR, 'teal/', i, ".teal.tellength.tsv")) }
  else if (type == "neo"){ as.character(paste0(I_DIR, 'neo/data_neo/', i, "/", i, ".neo.neoepitope.tsv.gz")) }
  else if (type == "neo_pep"){ as.character(paste0(I_DIR, 'neo/data_neo/', i,"/", i, ".neo.peptide_scores.tsv.gz")) }
}
reader <- function( i_file = "file_path", sample = "ACTN01020001T"){
    fread( i_file ) %>% mutate(sampleId = sample) 
}
