HELPER_DIR <- paste0(getwd(),'/helpers/')

source(paste0(HELPER_DIR,'map.r'))
source(paste0(HELPER_DIR,'shortcuts.r'))
source(paste0(HELPER_DIR,'helpers.r'))

go <- 
fread(paste0("data/0_run.csv")) %>% 
 fi(!grepl("rna_geneset", x), grepl("driver", x), !grepl("pathway", x), x != "drivers_total") %>%
 #fi(grepl("hotspot", x)) %>%
 fi(cohort != "Unknown primary (e.g. CUP)")

head(go)

options(repr.plot.width = 10, repr.plot.height = 6)

go %>%
 fi(y == "teal_final_ratio", x %in% c("driver_MECOM", "driver_CREBBP", "driver_NTHL1", "driver_ATRX", "driver_POT1", "driver_ARID1B"), covariate == "+ as.factor(biopsy) + purity + purity_ploidy") %>%
 ar(est) %>% 
 ggplot( aes(x = est, y = reorder(cohort, est))) + 
 geom_point() + 
 geom_errorbarh(aes(xmin = est - 2*se, xmax = est + 2*se), height = 0.2) +
 theme_bw() + 
 facet_wrap(~x, scales = "free_x", ncol = 6)

base <- fread(paste0(SHARE_DIR, "biomarkers_base.csv"))
base <- rbind(base, base %>% mu(cohort = "Pan-Cancer")) %>% gb(cohort) %>% mu(ct = n()) %>% fi(ct > 40, cohort != "Unknown primary (e.g. CUP)")

base_ready <- 
base %>% 
 se(sampleId, cohort, contains("teal"), any_of( unique(go %>% pu(x))))  %>%
 ga( x, val, -sampleId, -cohort, -teal_ref_raw,-teal_tumor_raw, -teal_tumor_final,-teal_raw_ratio,	-teal_final_ratio) %>%
 gb( cohort, x ) %>%
 mu( total_events = sum(val)) %>% 
 fi( total_events >= 3 )

fixer <- base %>% gb(cohort) %>% su(total_patients = n())

idx <- base_ready %>% se(cohort, x, total_events) %>% unique()

options(repr.plot.width = 10, repr.plot.height = 3)

teal_outputs_vs_ploidy <- 
base %>% 
 se( sampleId, cohort, purity_ploidy, contains("teal") ) %>%
 fi( cohort == "Pan-Cancer") %>% 
 ga( teal_summary, val, -purity_ploidy, -cohort, -sampleId) %>%
 fi( teal_summary != "teal_tumor_raw") %>% 
 ggplot( aes(x = purity_ploidy, y = val)) + 
 geom_point() + 
 facet_wrap(~teal_summary, scales = "free", ncol = 4) + 
 theme_bw() + 
 labs(title = "Telomere Length Measurement vs Ploidy", y = "Telomere Length Measurement", x = "Ploidy") + 
 theme(plot.title = element_text(hjust = .5)) + 
 geom_smooth()

teal_outputs_vs_ploidy

all_patients <- base %>% gb(cohort) %>% su(total_patients = n())

co_occurence_atrx <- 
base %>% 
 se( sampleId, cohort, contains("driver"), -contains("pathway"), -drivers_total ) %>%
 ga( x, val, -sampleId, -cohort, -driver_ATRX) %>%
 gb( cohort, x) %>%
 su( total_driver_events = sum(val), 
     correlation_ATRX = cor(driver_ATRX, val, use = "pairwise.complete.obs"), 
     total_ATRX = sum(driver_ATRX),
     total_cooccurrence_ATRX = sum(driver_ATRX + val == 2)) %>%
 ug()

co_occurence_tert <- 
base %>% 
 se( sampleId, cohort, contains("driver"), -contains("pathway"), -drivers_total ) %>%
 ga( x, val, -sampleId, -cohort, -driver_TERT) %>%
 gb( cohort, x) %>%
 su( total_driver_events = sum(val), 
     correlation_TERT = cor(driver_TERT, val, use = "pairwise.complete.obs"), 
     total_TERT = sum(driver_TERT),
     total_cooccurrence_TERT = sum(driver_TERT + val == 2)) %>%
 ug()

co_occurence_pot1 <- 
base %>% 
 se( sampleId, cohort, contains("driver"), -contains("pathway"), -drivers_total ) %>%
 ga( x, val, -sampleId, -cohort, -driver_POT1) %>%
 gb( cohort, x) %>%
 su( total_driver_events = sum(val), 
     correlation_POT1 = cor(driver_POT1, val, use = "pairwise.complete.obs"), 
     total_POT1 = sum(driver_POT1),
     total_cooccurrence_POT1 = sum(driver_POT1 + val == 2)) %>%
 ug()

cocurrence <- 
co_occurence_atrx %>%
 full_join(co_occurence_pot1 , by = c("cohort", "x", "total_driver_events")) %>%
 full_join(co_occurence_tert , by = c("cohort", "x", "total_driver_events")) %>%
 full_join(all_patients, by = "cohort") %>% 
 mu(x = gsub("driver_", "", x))

go_ready <- 
go %>% 
 ij(idx, by = c("cohort", "x")) %>%
 mu(x = gsub("driver_", "", x)) %>%
 lj(cocurrence, by = c("cohort", "x")) %>%
 se(-type, -data, -model, -lrt_pval) %>%
 tm(cohort, total_patients, measurement = y, driver = x, total_driver_events, covariate, est, se, pval, 
    correlation_ATRX = round(correlation_ATRX, 2), 
    correlation_POT1 = round(correlation_POT1,2), 
    correlation_TERT = round(correlation_TERT, 2), 
    total_ATRX, total_POT1, total_TERT,	
    total_cooccurrence_ATRX, total_cooccurrence_POT1, total_cooccurrence_TERT,
    frac_cooccurrence_ATRX = round(total_cooccurrence_ATRX/total_ATRX, 2), 
    frac_cooccurrence_POT1 = round(total_cooccurrence_POT1/total_POT1, 2), 
    frac_cooccurrence_TERT = round(total_cooccurrence_TERT/total_TERT, 2)) %>%
 fi(measurement != "teal_ref_raw") %>%
 ug() %>% 
 fi( (cohort == "Pan-Cancer" & covariate == "+ as.factor(cohort) + as.factor(biopsy) + purity + purity_ploidy") | 
     (cohort != "Pan-Cancer" & covariate == "+ as.factor(biopsy) + purity + purity_ploidy"))

share <- 
go_ready %>%
 gb(driver) %>%
 mu(min_p = min(pval)) %>%
 fi(min_p < .001) %>%
 ar( desc(total_patients), pval)

fwrite(share, "data/share_raw_output.csv")

options(repr.plot.width = 10)

library(cowplot)
library(ggrepel)

plotter <- function( y, title ) {
go_ready %>% 
 fi(measurement == y, cohort == "Pan-Cancer", pval < .001) %>% 
 ggplot( aes( x = est, y = log2(-log10(pval)))) +
 geom_text_repel(aes(label = driver )) + 
 theme_bw() + 
 labs( y = "Log2( -Log 10( p-value ))", 
       x = paste0("Driver association with ", title), 
       title = paste0(title, " vs Driver Signals")) +
 theme(plot.title = element_text(hjust = .5))
}

plt_names <- 
c("teal_final_ratio" = "Tumor vs Germline Telomere Ration (Final)", 
  "teal_raw_ratio" = "Tumor vs Germline Telomere Ration (Raw)",
  "teal_tumor_final" = "Tumor Telomere Length (Final)",
  "teal_tumor_raw" = "Tumor Telomere Length (Raw)")

plts <- list()
for( i in names(plt_names)){
 plts[[i]] <- plotter(i, plt_names[i])
}

options(repr.plot.width = 11, repr.plot.height = 8)

overall <- 
plot_grid(plts$teal_final_ratio, 
          plts$teal_raw_ratio,
          plts$teal_tumor_final, 
          plts$teal_tumor_raw, 
          labels = "AUTO",   # Optional: Adds labels A, B, C, D
          ncol = 2,          # Number of columns
          align = "hv")

overall

ggsave("data/overall.png", overall, width = 11, height = 8)

options(repr.plot.width = 15, repr.plot.height = 9)

cohort_level <-
go_ready %>% 
 fi(cohort != "Pan-Cancer", pval < .01, grepl("ratio", measurement)) %>% 
 ggplot( aes( x = est, y = -log10(pval), color = measurement)) +
 geom_text_repel( aes(label = driver ), size = 3) + 
 facet_wrap(~cohort, ncol = 7, scales = "free") + 
 theme_bw() + 
 scale_x_continuous(expand = expansion(mult = c(.4, .4))) +
 scale_y_continuous(expand = expansion(mult = c(.4, .4))) + 
 labs( y = "-Log 10( p-value )", 
       x = "Effect Size", 
       title = "Telomere Length Ratio (Tumor vs Germline) vs Driver Events)") + 
 theme(plot.title = element_text(hjust = .5))

cohort_level

ggsave("data/cohorts.png", cohort_level, width = 15, height = 9)
