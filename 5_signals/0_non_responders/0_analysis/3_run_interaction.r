source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "fisher.r"))

library(cluster)
library(purrr)

univariate_results <- fread(paste0(SHARE_DIR, "2_run_marginal_output.csv")) 
fisher_base <- fread(paste0(SHARE_DIR, "fisher_base.csv"))

treatment_mechanism_map <- fread(paste0(SHARE_DIR, "treatment_mechanism_map.csv"))

top <- univariate_results %>% fi(fisher_pval < .01, surv_pval < .01, direction == "Non-Response") 

auto_cluster_features <- function(data, method = "average", max_clusters = 5, corr_method = "pearson") {
  # Ensure features are columns
  if (nrow(data) < ncol(data)) {
    warning("Data has more features than samples; transposing for correlation.")
    data <- t(data)
  }
  
  # Step 1: Compute correlation matrix and convert to distance
  corr_matrix <- cor(data, method = corr_method, use = "pairwise.complete.obs")
  dist_matrix <- as.dist(1 - abs(corr_matrix))
  
  # Step 2: Hierarchical clustering
  hc <- hclust(dist_matrix, method = method)
  
  # Step 3: Evaluate silhouette scores for different cluster numbers
  best_score <- -Inf; 
  best_k <- 2; 
  best_clusters <- NULL

  for (k in 2:min(max_clusters, ncol(data) - 1)) {
    clusters <- cutree(hc, k = k)
    sil <- silhouette(clusters, dist_matrix)
    avg_sil <- mean(sil[, "sil_width"])
    
    if (avg_sil > best_score) {
      best_score <- avg_sil
      best_k <- k
      best_clusters <- clusters
    }
  }
  cat("Best number of clusters:", best_k, "\n")
  return(best_clusters)
}

add_cluster_labels <- function(i = "Skin Melanoma ## Immunotherapy"){
  data <- 
    fisher_base %>% 
      fi(cohortGo == i) %>% 
      se(any_of(top %>% fi(cohortGo == i) %>% pu(feature)))

  data.frame(auto_cluster_features(data)) %>% 
   rownames_to_column("feature") %>% 
   rename(cluster = auto_cluster_features.data.) %>% 
   mu(cohortGo = i)
}

cluster_labels <- data.frame()
for( i in unique(top$cohortGo)){
  print(i); flush.console();
  tmp <- tryCatch({add_cluster_labels(i)}, error = function(e) {return(NA)})
  if(is.data.frame(tmp)) cluster_labels <- rbind(cluster_labels, tmp)  
}

top_go <- 
top %>% 
 lj(cluster_labels, by = c("cohortGo", "feature")) %>% 
 gb(cohortGo, cluster) %>% mu(rk = row_number(fisher_pval)) %>% fi(rk <= 5) %>%
 se(cohortGo, feature, cluster) %>% 
 ug() %>% 
 drop_na()

features <- top_go %>% fi(i == cohortGo) %>% pu(feature)
clusters <- top_go %>% fi(i == cohortGo) %>% pu(cluster)

feature_pair <- list()
for( i in unique(top$cohortGo)){
    features <- top_go %>% fi(i == cohortGo) %>% pu(feature)
    clusters <- top_go %>% fi(i == cohortGo) %>% pu(cluster)
    if(length(features) > 1){
        feature_pair[[i]][["pairs"]] <- combn(features, 2, simplify = FALSE)
        feature_pair[[i]][["clusters"]] <- combn(clusters, 2, simplify = FALSE)
    } 
}

cluster_index <- data.frame()
for(i in names(feature_pair)){
  tmp_and <- 
  df(cohortGo = i,
     feature = unlist(lapply(feature_pair[[i]]$pairs, function(i) paste0(i[1], "_and_", i[2]))),
     clusters = unlist(lapply(feature_pair[[i]]$clusters, function(i) length(unique(i)))))

  tmp_or <- 
  df(cohortGo = i, 
     feature = unlist(lapply(feature_pair[[i]]$pairs, function(i) paste0(i[1], "_or_", i[2]))),
     clusters = unlist(lapply(feature_pair[[i]]$clusters, function(i) length(unique(i)))))
                      
  cluster_index <- rbind(cluster_index, tmp_and)
  cluster_index <- rbind(cluster_index, tmp_or)                    
}

fwrite(cluster_index, paste0(SHARE_DIR, "cluster_index.csv"))

add_combination_feature <- function(i = "Skin Melanoma ## Immunotherapy", pair = c('clin_hasRadiotherapyPreTreatment','driver_B2M')){

  fisher_base %>% 
   fi( cohortGo == i) %>% 
   se( sampleId, cohortGo, non_response, pfsEvent, daysToPfsEvent, primaryTumorLocation, purity, any_of(pair)) %>% 
   mu( 
    !!paste0(pair[1], "_and_", pair[2]) := 
     case_when(
      is.na(!!sym(pair[1])) | is.na(!!sym(pair[2])) ~ NA,
      (!!sym(pair[1]) + !!sym(pair[2])) == 2 ~ 1,
      TRUE ~ 0)#, 
#  !!paste0(pair[1], "_or_", pair[2]) := 
#     case_when(
#      is.na(!!sym(pair[1])) | is.na(!!sym(pair[2])) ~ NA,
#      (!!sym(pair[1]) + !!sym(pair[2])) >= 1 ~ 1,
#      TRUE ~ 0) 
  ) %>%  
  se(-any_of(pair))
}

cohort_combinations <- function(i = "Skin Melanoma ## Immunotherapy"){
  pairs <- feature_pair[[i]][['pairs']]
    
  combos <- list()
  for( j in seq(length(pairs)) ){
    pair_ready <- pairs[[j]]
    combos[[j]] <- add_combination_feature(i = i, pair = pair_ready)
  }
  ## Do separately for same or different clusters   
  clean_combos <- Filter(Negate(is.null), combos)  
  reduce(clean_combos, ~ inner_join(.x, .y, by = c("sampleId", "cohortGo", "non_response", "pfsEvent", "daysToPfsEvent", "primaryTumorLocation", "purity")))  
} 

combos_ready <- list()

for( i in names(feature_pair)){
    print("Finding interactions: "); print(i); flush.console()
    if( !i %in% names(combos_ready)) {
      combos_ready[[i]] <- cohort_combinations(i = i)
    }
}

saveRDS(combos_ready, paste0(TMP_DIR, "combo_features.Rds"))

fisher_combos <- list()
for( i in names(combos_ready)){
  print(i); flush.console(); 
  fisher_combos[[i]] <- tryCatch({ra_formatter_and_test(combos_ready[[i]] %>% se(-sampleId, -pfsEvent, -daysToPfsEvent, -primaryTumorLocation, -purity))}, error = function(e){NA})
}

scanner <- function (y = "Surv(daysToPfsEvent, pfsEvent)", features, covariates, df = "df", mod = "coxph", cohort = "Pan-Cancer") {
    out <- data.frame()
    for (f in features) {
        if (grepl("rna_", f)) {tmp <- get_stats(y = y, x = f, covariate = paste0(covariates, "+ purity"), data = df, model = mod)} 
        else {tmp <- get_stats(y = y, x = f, covariate = covariates, data = df, model = mod)}
        if (is.data.frame(tmp)) out <- rbind(out, tmp)
    }
    out
}

pfs_results <- list()
for( i in names(combos_ready)){
    print(i); flush.console(); 
    tmp <- combos_ready[[i]] %>% se(-sampleId)
    categorical_features <- names(tmp %>% se(-cohortGo, -non_response, -pfsEvent, -daysToPfsEvent, -primaryTumorLocation, -purity))
    pfs_results[[i]] <- scanner(  feature = categorical_features, covariates = "", df = "tmp") %>% mu(cohortGo = i)
}

combo_results <- data.frame()
for( i in names(combos_ready)){
    print(i); flush.console();
    fisher_i <- fisher_combos[[i]]
    pfs_i <-  pfs_results[[i]]
    add <- 
    fisher_i %>% 
     lj(pfs_i %>% tm(cohortGo, feature = x, surv_est = est, surv_se = se, surv_pval = pval), 
        by = c("cohortGo", "feature")) 
    combo_results <- rbind(combo_results, add)
}

combo_results_ready <- 
combo_results %>% 
 mu(type = "combination", dcb = responders, no_dcb = non_responders, ct = no_dcb + dcb) %>% 
 lj(univariate_results %>% se(cohortGo, group) %>% unique(), by = "cohortGo")

nr_threshold <- .05
fdr_threshold <- .05

combiner <- function(a){
 b <- strsplit(a, "_")[[1]]
 paste0(b[-length(b)], collapse = "_")
}

together <- 
rbind(univariate_results %>% mu(type = "univariate"), combo_results_ready) %>% 
 lj(cluster_index, by = c("cohortGo", "feature")) %>% 
 mu(clusters = ifelse(is.na(clusters), 1, clusters)) %>% 
 rw() %>% mu(short_feature = combiner(feature)) %>% ug() %>% 
 gb(cohortGo, short_feature, fisher_pval) %>% mu(tot = n()) %>% ug() %>% 
 #fi(group != "type") %>%  
 fi(tot == 1 | direction == "Non-Response") %>% 
 mu(or = ifelse(or == "Inf", exp(5), or), p_bh = p.adjust(fisher_pval, method = "fdr")) %>% 
 rw() %>% mu( derived_treatmentName = str_split_fixed( cohortGo	, " ## ", n = 2)[2]) %>% ug() %>% 
 lj( treatment_mechanism_map , by = "derived_treatmentName") %>% 
 mu( derived_treatmentMechanism = ifelse(is.na(derived_treatmentMechanism), derived_treatmentName, derived_treatmentMechanism),
     treatment = derived_treatmentName, mechanism = gsub(" ## ", "/", derived_treatmentMechanism),
     alpha_gp = case_when(p_bh >= fdr_threshold ~ "none", p_bh < nr_threshold ~ "dark"), 
     highlight = ((e_nr/events >= 1 - nr_threshold) & 
                  (p_bh < fdr_threshold) & 
                  !((type == "combination") & (clusters == 1))),
     super_highlight = highlight & surv_pval < .01, 
     color_gp = case_when(highlight ~ mechanism, !highlight & (p_bh < fdr_threshold) ~ "significant",TRUE ~ "none")) %>% 
 fi(!grepl("battle", feature), treatment != "Immunotherapy")

fwrite(together %>% fi(super_highlight), paste0(SHARE_DIR,"share_highlights_with_fran_update.csv"))

fwrite( together, paste0(SHARE_DIR, "3_run_interaction_combined_output.csv"))
