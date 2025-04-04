source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

somatic_exome <- 
fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/somatic_exome.csv") %>% 
 gb(chromosome, position, gene, tier) %>% 
 mu(ct = n(), hotspot = (tier == "HOTSPOT" | ct > 30)) %>% 
 ug()

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

#som_drivers <- rbind(hotspots, biallelic, inframe) %>% lj(cn_gene %>% se(gene, chromosomeBand), by = "gene")

sel_cvs <- fread("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/dnds/sel_cvs.csv")

w_thresh <- 2
q_thresh <- .05

gper <- function(i){
 if(grepl("non", i)){ "nonsense" }
 else if(grepl("mis", i)){ "missense" }       
 else { "indel" }
}

dnds_drivers_ref <- 
sel_cvs %>% 
 fi(grepl("Pan-Cancer", cohort), 
    (qtrunc_cv < q_thresh | qmis_cv < q_thresh | qind_cv < q_thresh),
    (wmis_cv > w_thresh | wnon_cv > w_thresh | wind_cv > w_thresh )) %>% 
 mu(tot = n_syn	+ n_mis	+ n_non + n_spl	+ n_ind) %>% 
 se(gene_name, tot, wmis_cv, wnon_cv, wind_cv, qmis_cv, qtrunc_cv, qind_cv) %>% 
 ga( w, wval, -gene_name, -tot, -qmis_cv, -qtrunc_cv, -qind_cv) %>% 
 ga( q, qval, -gene_name, -tot, -w, -wval) %>% 
 fi(wval > w_thresh, qval < q_thresh, 
   (grepl("mis", w) & grepl("mis", q)) | (grepl("ind", w) & grepl("ind", q)) | (grepl("non", w) & grepl("trunc", q))) %>%
 rw() %>% mu(type = gper(w)) %>% ug() %>% 
 tm(gene = gene_name, type)

annotater <- function(i){
 if(grepl("stop", i)){ "nonsense" }
 else if(grepl("missense", i)){ "missense" }       
 else if(grepl("inframe", i)){ "indel" }
 else if(grepl("frameshift", i)){ "indel" }    
 else if(grepl("start_lost", i)){ "start_lost" }
 else {"synonymous"}
}

dnds <- 
somatic_exome %>% 
 fi( gene %in% unique(dnds_drivers_ref$gene)) %>% 
 rw() %>% 
 mu( type= annotater(annotation)) %>% 
 ug() %>% 
 fi( type != "synonymous") %>% 
 se( sampleId, chromosome, gene, transcript, type, biallelic ) %>% 
 ij( dnds_drivers_ref, by = c("gene", "type")) %>% 
 mu( method = "DNDS" ) %>% 
 rename(annotation = type)

somatic_drivers <- rbind(dnds, hotspots, biallelic, inframe) 

fwrite(somatic_drivers, "/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/database/drivers_full/somatic.txt")
