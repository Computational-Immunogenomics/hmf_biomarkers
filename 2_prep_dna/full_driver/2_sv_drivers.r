source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

som_exome <- fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/somatic_exome.csv")

breakends <- fread( paste0(TMP_DIR, "structural_variants/breakend.csv")) 

drivers <- fread( paste0(TMP_DIR, "drivers.csv"))

amps <- fread( "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/drivers_full/amp.txt")

dels <- fread( "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/drivers_full/del.txt")

dels_cts <- dels %>% gb(gene) %>% su(ct = n()) %>% ar(desc(ct)) %>% fi(ct < 1000)

amps_cts <- amps %>% gb(gene) %>% su(ct = n()) %>% ar(desc(ct)) %>% fi(ct < 1000)

amps_cts %>% lj(drivers %>% tm(gene, present = 1) %>% unique(), by = "gene")

cts %>% fi(grepl("IFN", gene))

hist(cts$ct)

amps %>% gb(gene) %>% su(ct = n()) %>% ar(desc(ct))

disruptions <- 
breakends %>%
 fi(disruptive) %>% 
 tm(sampleId, chromosome, chromosomeBand = chrBand, gene, transcriptId, canonical, disruptive, method = "DISRUPTION") %>% 
 unique()

head(disruptions)

disruptions %>% 
 gb(chromosome, chromosomeBand, gene) %>% 
 su(ct = n()) %>% 
 ar(desc(ct))
