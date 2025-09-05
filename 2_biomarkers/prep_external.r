source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

smoking <- fread("/mnt/bioinfnas2/immunocomp/manuel/tme/0_organise_data/output_data/aetiology_features/smoker_status.csv.gz")
colibactin <- fread("/mnt/bioinfnas2/immunocomp/manuel/tme/0_organise_data/output_data/aetiology_features/colibactin_status.csv.gz")
signatures <- fread("/mnt/petasan_immunocomp/projects/CUPs/processed_data/mut_signatures/signatures_all_cancers.csv.gz")

smoking_ready <- smoking %>% tm(sampleId, clin_smoker = smoker)

colibactin_ready <- colibactin %>% rename_with(~ paste0("colibactin_", .), -1)

sigs_ready <- 
signatures %>%
 se(-tissue_type) %>% 
 rename_with(~ paste0("signature_log_", .), -1) %>% 
 mu(across(.cols = -all_of("sampleId"), .fns = log1p))

external_ready <- 
smoking_ready %>% 
 full_join(colibactin_ready, by = "sampleId") %>% 
 full_join(sigs_ready, by = "sampleId")

fwrite(external_ready, paste0(READY_DIR, "external_ready.csv"))
