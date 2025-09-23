source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

print("Dimension of all biomarkers derived")
all <- fread(paste0(SHARE_DIR, "biomarkers_base.csv")); dim(all)

print("Dimension of biomarkers with purity, outcomes, removed features")
ready <- 
all %>% 
 fi(!is.na(purity), !is.na(durableClinicalBenefit)) %>% 
 mu(non_response = abs(durableClinicalBenefit-1)) %>%
 se(-contains("sv_")); dim(ready)

fwrite( ready %>% se(derived_treatmentName, derived_treatmentMechanism) %>% unique(), paste0(SHARE_DIR, "treatment_mechanism_map.csv"))

print("Dimension of RNA Data available")
dim(ready %>% fi(!is.na(rna_geneset_HALLMARK_MTORC1_SIGNALING)))

base_features <- 
names(
 ready %>% 
    se( contains("cider_"), contains("clin_"), contains("cn_"), contains("driver"), 
        contains("fusion"), contains("gie_"), contains("lilac_"), contains("neo_"), 
        contains("chord"), contains("purity"), contains("rna_"), contains("teal_"), 
        contains("viral_"), contains("hotspot"), contains("bacterial"), contains("teal"), 
        contains("signature")))

print("Number of initial features in the analysis")
length(base_features)

cn_features <- names(ready %>% se(contains("cn_"), -contains("lilac_")))

cn_ploidy_features <- ready %>% mutate(across(cn_features, ~ .x / purity_ploidy)) %>% se(all_of(cn_features)) 

cn_dels <- cn_ploidy_features %>% mu(across(everything(), ~ (. <= .7))) %>% rename_with(~ paste0(.x, "_del"))
cn_no_dels <- cn_ploidy_features %>% mu(across(everything(), ~ (. > .7))) %>% rename_with(~ paste0(.x, "_no_del"))
cn_amps <- cn_ploidy_features %>% mu(across(everything(), ~ (. >= 1.3))) %>% rename_with(~ paste0(.x, "_amp"))
cn_no_amps <- cn_ploidy_features %>% mu(across(everything(), ~ (. < 1.3))) %>% rename_with(~ paste0(.x, "_no_amp"))

tmb_msi <- 
data.frame(
    exp(ready %>% se(contains("PerMb"))-1), 
        ready %>% se(purity_msStatus, purity_tmbStatus, purity_tmlStatus, chord_hrStatus))

tmb_features <- 
tmb_msi %>% 
 mu( purity_tmbPerMb_gt2 = (purity_tmbPerMb > 2),
     purity_tmbPerMb_gt4 = (purity_tmbPerMb > 4),
     purity_tmbPerMb_gt6 = (purity_tmbPerMb > 6),
     purity_tmbPerMb_gt8 = (purity_tmbPerMb > 8),
     purity_tmbPerMb_gt12 = (purity_tmbPerMb > 12), 
     purity_tmbPerMb_lt2 = (purity_tmbPerMb < 2),
     purity_tmbPerMb_lt4 = (purity_tmbPerMb < 4),
     purity_tmbPerMb_lt6 = (purity_tmbPerMb < 6),
     purity_tmbPerMb_lt8 = (purity_tmbPerMb < 8),
     purity_tmbPerMb_lt12 = (purity_tmbPerMb < 12), 
     purity_tmlStatus_high = purity_tmlStatus, 
     purity_tmlStatus_low = 1 - purity_tmlStatus, 
     purity_tmbStatus_high = purity_tmbStatus, 
     purity_tmbStatus_low = 1 - purity_tmbStatus, 
     chord_hrStatus_high = chord_hrStatus, 
     chord_hrStatus_low = 1 - chord_hrStatus
   ) %>% 
 se(contains("purity_tmbPerMb_"), contains("Status_"))

msi_features <- 
tmb_msi %>% 
  mu(purity_msIndelsPerMb_gt1 = (purity_msIndelsPerMb > 1),
     purity_msIndelsPerMb_gt2 = (purity_msIndelsPerMb > 2),
     purity_msIndelsPerMb_lt1 = (purity_msIndelsPerMb < 1),
     purity_msIndelsPerMb_lt2 = (purity_msIndelsPerMb < 2),
     purity_msStatus_high = purity_msStatus, 
     purity_msStatus_low = 1 - purity_msStatus, 
   ) %>% 
 se(contains("purity_msIndelsPerMb_"), contains("Status_"))

pathways_affected <- 
 ready %>% 
    mu( driver_pathways_lt3 = as.numeric(drivers_pathway_total < 3), 
        driver_pathways_lt5 = as.numeric(drivers_pathway_total < 5),
        driver_pathways_lt7 = as.numeric(drivers_pathway_total < 7),
        driver_pathways_lt9 = as.numeric(drivers_pathway_total < 9), 
        driver_pathways_gt1 = as.numeric(drivers_pathway_total > 1), 
        driver_pathways_gt2 = as.numeric(drivers_pathway_total > 2), 
        driver_pathways_gt4 = as.numeric(drivers_pathway_total > 4),
        driver_pathways_gt6 = as.numeric(drivers_pathway_total > 6),
        driver_pathways_gt8 = as.numeric(drivers_pathway_total > 8)) %>%  
  se(contains("driver_pathways_lt"), contains("driver_pathways_gt"))

rest <- 
 ready %>% 
   se(-all_of(cn_features), 
      -purity_msStatus, 
      -purity_tmbStatus, 
      -purity_tmlStatus, 
      -chord_hrStatus, 
      -drivers_pathway_total, 
      -contains("PerMb") )

filter_ref <- 
data.frame(
   "mn" = apply(rest %>% se(any_of(base_features)), 2, mean, na.rm = TRUE), 
   "sd" = apply(rest %>% se(any_of(base_features)), 2, sd, na.rm = TRUE), 
   "zeros" = apply(rest %>% se(any_of(base_features)) == 0, 2, sum, na.rm = TRUE),
   "non_zeros" = apply(rest %>% se(any_of(base_features)) != 0, 2, sum, na.rm = TRUE), 
   "nas" = apply(is.na(rest %>% se(any_of(base_features))), 2, sum, na.rm = TRUE), 
   "non_nas" = apply(!is.na(rest %>% se(any_of(base_features))), 2, sum, na.rm = TRUE)) %>% 
 mu(pct_zeros = zeros/(zeros+non_zeros), pct_nas = nas/(nas+non_nas)) %>%
 rownames_to_column("feature") 

binary_features <- 
rest %>% 
 se(all_of(filter_ref %>% pu(feature))) %>% 
 se(where(~all(. %in% c(0, 1, NA))))

non_binary_non_sparse_features <- 
rest %>% 
 se(all_of(filter_ref %>% fi(pct_zeros <= .5) %>% pu(feature))) %>% 
 se(!where(~all(. %in% c(0, 1, NA)))) 

non_binary_sparse_features <- 
rest %>% 
 se(all_of(filter_ref %>% fi(pct_zeros > .5) %>% pu(feature))) %>% 
 se(!where(~all(. %in% c(0, 1, NA)))) 

integerer <- function(df){
 df[] <- lapply(df, function(x) if(is.logical(x)) as.integer(x) else x); df    
}

binarify <- function(df, threshold = 50, direction = "gt" ){
 if(direction == "gt"){
   tmp <- df %>% 
    mu(across(everything(), ~ (. >= quantile(., threshold/100, na.rm = TRUE)))) %>% 
    rename_with(~ paste0(.x, "_gt", as.character(threshold)))
 } else {
   tmp <- df %>% 
    mu(across(everything(), ~ (. < quantile(., threshold/100, na.rm = TRUE)))) %>% 
    rename_with(~ paste0(.x, "_lt", (as.character(threshold))))
 }
 integerer(tmp)
}

binarify_sparse <- function(df, direction = "gt" ){
 if(direction == "gt"){
  tmp <- df %>% mu(across(everything(), ~(. > 0))) %>% rename_with(~ paste0(.x, "_gt0"))
 } else {
  tmp <- df %>% mu(across(everything(), ~ (. == 0))) %>% rename_with(~ paste0(.x, "_eq0"))
 }
 integerer(tmp)
}

gt50 <- binarify(non_binary_non_sparse_features, 50, "gt")
lt50 <- binarify(non_binary_non_sparse_features, 50, "lt")

smooth_features <- 
unlist(lapply(names(Filter(function(x) .48 < x && x < .52, apply(gt50, 2, mean, na.rm = TRUE))), 
                           function(i) strsplit(i, "_gt50")[[1]][1]))

gt25 <- binarify(non_binary_non_sparse_features %>% se(any_of(smooth_features)), 25, "gt")
gt75 <- binarify(non_binary_non_sparse_features %>% se(any_of(smooth_features)), 75, "gt")
lt25 <- binarify(non_binary_non_sparse_features %>% se(any_of(smooth_features)), 25, "lt")
lt75 <- binarify(non_binary_non_sparse_features %>% se(any_of(smooth_features)), 75, "lt")

sparse_gt <- binarify_sparse(non_binary_sparse_features, "gt")
sparse_lt <- binarify_sparse(non_binary_sparse_features, "lt")

prepared_features <- 
cbind(binary_features, 
      gt50, 
      lt50, 
      gt25, 
      gt75, 
      lt25, 
      lt75, 
      sparse_gt, 
      sparse_lt, 
      pathways_affected, 
      cn_dels, 
      cn_amps, 
      tmb_features, 
      msi_features) %>%
  mutate(across(where(is.logical), as.integer)) %>% 
  select(where(~is.numeric(.) && sum(., na.rm = TRUE) >= 5))

go <- 
cbind( ready %>% se(-contains(base_features), purity),  prepared_features) %>% 
  rw() %>%
  mu(groupedTreatmentType = paste0(unique(strsplit(derived_treatmentType, " ## ")[[1]]), 
     collapse = " ## ")) %>% 
  ug()

saveRDS( list("ready" = go, "features" = names(prepared_features)), paste0(SHARE_DIR, "biomarkers_ready.rds"))

print("Total Base Features")
print(length(base_features))

print("Base Features Counts")
table(unlist(lapply(base_features, function(i) strsplit(i, "_")[[1]][1])))

print("Derived Features counts")
table(unlist(lapply(names(prepared_features), function(i) strsplit(i, "_")[[1]][1])))

print("Total Number of Features")
print(length(names(prepared_features)))

print("Total RNA Features")
rna_features <- sum(unlist(lapply(names(prepared_features), function(i) strsplit(i, "_")[[1]][1])) == "rna")
rna_features

print("Total Non-RNA Features")
print(length(names(prepared_features)) - rna_features)
