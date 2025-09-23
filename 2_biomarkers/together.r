source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

setwd(READY_DIR)

sources <- 
c("purity",
  "clinical",
  "drivers", 
  "cn_simple",
  "isofox_genesets",
  "cider_dna", 
  "teal", 
  "viral", 
  "lilac", 
  "neo", 
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
  dfs[[i]] <- fread( paste0(i, "_ready.csv")) %>% mu(sampleId = as.character(sampleId))
})

ready <- 
 dfs %>% 
 reduce(left_join, by = "sampleId") %>%
 rename_with(~ gsub("[^A-Za-z0-9_]", "_", .x))

ready <- ready %>% mutate( across(c(contains("driver"), contains("viral"), contains("hotspot")), ~replace_na(.,0) ))

fwrite( ready, paste0(SHARE_DIR, "biomarkers_base.csv") )
