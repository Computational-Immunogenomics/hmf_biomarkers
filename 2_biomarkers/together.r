source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

setwd(READY_DIR)

sources <- 
c("clinical", 
  "purity", 
  "drivers", 
  "cn_simple",
  "isofox_genesets",
  "cider_dna", 
  "teal", 
  "viral", 
  "lilac", 
  "neo", 
  "neo_pep", 
  "svs", 
  "gie", 
  "fusions_dna", 
  "chord",
  "hotspots", 
  "external")

dfs <- list()
system.time(
for( i in sources){
  print(i); flush.console()
  dfs[[i]] <- fread( paste0(i, "_ready.csv"))
})

for( i in names(df)) {
    print(i); print(dim(dfs[[i]]))
}

ready <- 
 dfs %>% 
 reduce(left_join, by = "sampleId") %>%
 rename_with(~ gsub("[^A-Za-z0-9_]", "_", .x)) %>% 
 se(-contains("perplexity"))

ready <- 
ready %>% mutate( across(c(contains("driver"), contains("viral"), contains("hotspot")), ~replace_na(.,0) ))

fwrite( ready, paste0(SHARE_DIR, "biomarkers_base.csv") )
