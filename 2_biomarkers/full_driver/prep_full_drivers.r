source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

som_exome <- fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/somatic_exome.csv")

drivers <- fread( paste0(TMP_DIR, "drivers.csv"))

#dnds <- readRDS("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/dnds/dnds.Rds")
sel_cvs <- fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/dnds/sel_cvs.csv")

sig <- 
sel_cvs %>% 
 fi(grepl("Pan-Cancer", cohort)) %>% 
 fi(qglobal_cv < .2 | qtrunc_cv < .2 | qmis_cv < .2 | qind_cv < .2) %>% 
 fi(wmis_cv > 2 | wnon_cv > 2 | wind_cv > 2 | qtrunc_cv > 2 ) %>% 
 #fi( gene_name %in% genes) %>% 
 #fi(gene_name %in% c(drivers %>% pu(gene) %>% unique())) %>% 
 mu(tot = n_syn	+ n_mis	+ n_non + n_spl	+ n_ind) %>% 
 fi(tot > 10) %>% 
 se(gene_name, tot, wmis_cv, wnon_cv, wind_cv, qmis_cv, qtrunc_cv, qind_cv, qglobal_cv, pglobal_cv, cohort) %>% ar(qglobal_cv)

sig %>% fi(grepl("TGFB", gene_name))

purities <- fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/purities.csv")

cn_gene <- fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/cn_gene.csv")

breakends <- fread( paste0(TMP_DIR, "structural_variants/breakend.csv")) 

drivers <- fread( paste0(TMP_DIR, "drivers.csv"))

drivers %>% 
 fi(grepl("TGFB", gene)) %>% 
 gb(gene, likelihoodMethod) %>% 
 su(ct = n()) %>% 
 ar(desc(ct))

annotation <- fread(paste0(dirname(dirname(getwd())), "/util/annotation/pathways.csv"))

drivers <- fread( paste0(TMP_DIR, "drivers.csv"))

#dnds <- readRDS("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/dnds/dnds.Rds")
sel_cvs <- fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/dnds/sel_cvs.csv")

drivers %>% 
 tm(gene_name = gene, category) %>% 
 unique() %>% 
 rj( sel_cvs %>% fi(cohort == "Pan-Cancer") %>% se(-contains("p"), -cohort), by = "gene_name") %>%
 ar(qglobal_cv) %>% 
 mu(nm_ratio = wnon_cv/wmis_cv) %>% 
 mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
 relocate(nm_ratio, .after = category)

tmp <- 
sel_cvs %>% 
 mu(ct = n_syn + n_mis + n_non + n_spl + n_ind, 
    ratio = wnon_cv/wmis_cv, 
    type = ifelse(ratio > 1, "TSG", "ONC")) %>% 
 fi(cohort == "Pan-Cancer") %>% 
 se(cohort, gene_name, ct, type, ratio, contains("n_"), contains("w"), pglobal_cv, qglobal_cv) %>% 
 ar(desc(ct))

tmp

table(annotation$pathway)

tmp %>% 
 ar(qglobal_cv) %>% 
 fi(cohort ==  "Pan-Cancer") %>% 
 lj(drivers %>% tm(gene_name = gene, catalogue = TRUE, category) %>% unique(), by = "gene_name") %>% 
 se(catalogue, ct, gene_name, category, type, contains("w"), ratio, pglobal_cv, qglobal_cv) %>% 
 lj(annotation %>% rename(gene_name = gene), by = "gene_name") %>% 
 mu(rk = row_number()) %>% 
 fi(grepl("EGF", gene_name))

pan <- tmp %>% fi(cohort == "Pan-Cancer") %>% pu(gene_name) %>% unique()

sel_cvs %>% 
 fi(gene_name == "BMP5") %>% ar(pglobal_cv) %>% se(cohort, gene_name, contains("n_"), contains("w"), pglobal_cv) 

tmp %>% fi(!gene_name %in% pan, ct > 10) %>% ar(pglobal_cv)

dnds$globaldnds

see <- fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/drivers.csv")

drivers <- 
fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/drivers.csv") %>% 
 fi(likelihoodMethod == "DNDS") %>% 
 gb(gene) %>% 
 su(like = mean(driverLikelihood > .8), ct = n()) %>% mu(gene_name = gene) %>% 
 ug()

see %>% 
 gb(gene, likelihoodMethod) %>% 
 su(ct = n()) %>% 
 ar(desc(ct)) %>% 
 fi(gene == "PIK3CA")

dnds$sel_cv %>% ij(drivers, by = "gene_name") %>% se(gene_name, contains("n_"), contains("w"), like, ct, pglobal_cv) %>% gb(pglobal_cv < .01) %>% su(ct = n())

dnds$sel_cv %>% fi(qglobal_cv < .25) %>% fi(n_syn != 0)

somatic_exome <- 
som_exome %>% 
 gb(chromosome, position, gene, tier) %>% 
 mu(ct = n(), hotspot = (tier == "HOTSPOT" | ct > 30))

hotspots <- 
somatic_exome %>% 
 fi(hotspot) %>% 
 tm(sampleId, chromosome, gene, transcript, annotation, biallelic, method = "HOTSPOT")

biallelic <- 
somatic_exome %>% 
 fi(!hotspot, biallelic) %>% 
 tm(sampleId, chromosome, gene, transcript, annotation, biallelic, method = "BIALLELIC")

inframe <- 
somatic_exome %>% 
 fi(!hotspot, !biallelic, annotation %in% c("inframe_deletion", "inframe_insertion", "phased_inframe_insertion", "phased_inframe_deletion")) %>% 
 tm(sampleId, chromosome, gene, transcript, annotation, biallelic, method = "INFRAME")

som_drivers <- rbind(hotspots, biallelic, inframe) %>% lj(cn_gene %>% se(gene, chromosomeBand), by = "gene")

dim(cn_gene)

amps <- 
cn_gene %>% 
 ij(purities %>% se(sampleId, ploidy), by = "sampleId") %>%
 fi(maxCopyNumber > (3*ploidy)) %>% 
 mu( method = ifelse(minCopyNumber > 3*ploidy, "AMP", "PARTIAL_AMP"))

dels <- 
cn_gene %>% 
 ij( purities %>% se(sampleId, ploidy), by = "sampleId" ) %>%
 fi( minCopyNumber < .5 ) %>% 
 mu( method = "DEL" )  

length(unique(cn_gene$sampleId))

head(dels)

drivers %>% fi(grepl("TRAF3", gene))

dels %>% gb(chromosome, chromosomeBand, gene) %>% su(ct = n()) %>% ar(desc(ct)) %>% fi(chromosome == 23)

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
