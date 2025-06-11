source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "fisher.r"))

library(cluster)
library(purrr)

univariate_results <- fread(paste0(SHARE_DIR, "1_run_fishers_exact.csv")) 
fisher_base <- fread(paste0(SHARE_DIR, "fisher_base.csv"))

top <- univariate_results %>% fi(fisher_pval < combination_threshold) 

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
 gb(cohortGo, cluster) %>% mu(rk = row_number(fisher_pval)) %>% fi(rk <= 10) %>%
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

add_combination_feature <- function(i = "Skin Melanoma ## Immunotherapy", 
                                    pair = c('clin_hasRadiotherapyPreTreatment','driver_B2M')){
  fisher_base %>% 
   filter( cohortGo == i) %>% 
   select( sampleId, cohortGo, non_response, any_of(pair)) %>% 
   mutate( 
    !!paste0(pair[1], "_and_", pair[2]) := 
     case_when(
      is.na(!!sym(pair[1])) | is.na(!!sym(pair[2])) ~ NA,
      (!!sym(pair[1]) + !!sym(pair[2])) == 2 ~ 1,
      TRUE ~ 0), 
  !!paste0(pair[1], "_or_", pair[2]) := 
     case_when(
      is.na(!!sym(pair[1])) | is.na(!!sym(pair[2])) ~ NA,
      (!!sym(pair[1]) + !!sym(pair[2])) >= 1 ~ 1,
      TRUE ~ 0) 
  ) %>% se(-any_of(pair))
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
  reduce(clean_combos, ~ inner_join(.x, .y, by = c("sampleId", "cohortGo", "non_response")))  
} 

combos_ready <- list()

for( i in names(feature_pair)){
    print("Finding interactions: "); print(i); flush.console()
    if( !i %in% names(combos_ready)) {
      combos_ready[[i]] <- cohort_combinations(i = i)
    }
}

combo_results <- data.frame()

for( i in names(combos_ready)){
    print(i); flush.console()
    add <- tryCatch({ra_formatter_and_test( combos_ready[[i]] %>% se(-sampleId))}, error = function(e){NA})
    if(is.data.frame(add)){
        combo_results <- rbind(combo_results, add)
    }
}

combo_results <- data.frame()

for( i in names(combos_ready)){
    print(i); flush.console()
    add <- tryCatch({ra_formatter_and_test( combos_ready[[i]] %>% se(-sampleId))}, error = function(e){NA})
    if(is.data.frame(add)){
        combo_results <- rbind(combo_results, add)
    }
}

combo_results_ready <- 
combo_results %>% 
 mu(type = "combination", dcb = responders, no_dcb = non_responders, ct = no_dcb + dcb) %>% 
 lj(univariate_results %>% se(cohortGo, group) %>% unique(), by = "cohortGo")

together <- 
rbind(univariate_results %>% mu(type = "univariate"), combo_results_ready) %>% 
 lj(cluster_index, by = c("cohortGo", "feature")) %>% 
 mu(clusters = ifelse(is.na(clusters), 1, clusters))

fwrite( together, paste0(SHARE_DIR, "2b_go.csv"))
