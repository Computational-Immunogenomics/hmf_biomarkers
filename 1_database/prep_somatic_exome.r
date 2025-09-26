source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

EXOME_DIR <- paste0(TMP_DIR, "somatic_exome_big/")

remove <- c("intron_variant")

go <- data.frame()

system.time(
for( i in list.files(EXOME_DIR)){
  fp <- paste0(EXOME_DIR, i)
  print(fp); flush.console()
  tmp <- fread( fp ) %>% fi(! annotation %in% remove)
  go <- rbind(go, tmp)
})

fwrite(go, paste0(TMP_DIR, "somatic_exome.csv"))
