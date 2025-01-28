source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

rna_fusions <- fread(paste0(TMP_DIR, "isofox_fusions.csv"))

dna_fusions <- fread(paste0(TMP_DIR, "structural_variants/fusion.csv"))
dna_fusions2 <- fread(paste0(TMP_DIR, "structural_variants/vis_fusion.csv"))
