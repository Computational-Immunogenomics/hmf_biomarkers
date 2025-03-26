source("../hub/map.r")
source("../hub/shortcuts.r")
source("../hub/run.r")

#oo <- fread("1_run_output.csv")
oo <- fread("1_run_non_response_output.csv")

#head(oo)

base <- fread(paste0(O_DIR,"/biomarkers_go.csv")) %>% fi(!is.na(bestOverallResponse), !is.na(purity_purity))
ready <- readRDS("biomarker_cohorts2.Rds")

#head(oo)

names <- 
c("durableClinicalBenefit" = "Durable Benefit",
  "bestOverallResponse" = "Best Response",
  "Surv(daysToOsEvent, osEvent)" = "Overall Survival",
  "Surv(daysToPfsEvent, pfsEvent)" = "PFS Survival"
 )

oo <- 
oo %>% 
 mu(cohort = paste0( mechanism, "_", tissue, "_", type),
    p_adj_by = p.adjust(pval, method = "BY"), 
    #p_adj_bonf = p.adjust(pval, method = "bonferroni"), 
    p_adj_fdr = p.adjust(pval, method = "fdr"), 
    z = est/se) %>% 
 rw() %>% mu( y = names[[y]]) %>% ug()

low_event_counts <- function( type, mechanism, tissue = "all" ){
 
 if( tissue == "all"){  df <- ready[[type]][[mechanism]] } 
 else( df <- ready[[type]][[mechanism]][[tissue]] )
    
 df %>% 
  se(contains("driver_"), contains("viral_"), contains("hla_")) %>% 
  ga(x, event) %>% 
  mu(event = (event > 0)) %>% 
  gb(x) %>% 
  su(events = sum(event)) %>% 
  ug() %>% 
  fi(events < 10) %>% 
  mu(type = type, mechanism = mechanism, tissue = tissue) %>% 
  ar(desc(events))
}

remove <- data.frame()
for( i in c("exact", "contains", "tissue_exact", "tissue_contains")){
 for( j in names(ready[[i]])){ 
   if( i %in% c("exact", "contains")){
    tmp <- low_event_counts(i, j, "all")
    remove <- rbind(remove, tmp)
  } else {
    for( k in names(ready[[i]][[j]])){
       tmp <- low_event_counts(i, j, k)
       remove <- rbind(remove, tmp) 
    }   
   }
 } 
}

oo_go <- 
oo %>% 
 lj(remove, by = c("x", "type", "mechanism", "tissue")) %>% 
 fi(is.na(events), abs(z) < 10, abs(est) < 6) %>% 
 ar(pval) %>% 
 fi(mechanism != "Anti-PD") %>% 
 mu(tissue = ifelse(tissue == "all", "All", tissue), 
    type = ifelse(grepl("contains", type), "contains", "exact"),
    cohort2 = paste0(mechanism, "\n", tissue))

cts <- oo_go %>% gb(cohort2) %>% su(ct2 = max(n_samples))

oo_go <- 
oo_go %>% 
 lj(cts, by = "cohort2") %>% 
 mu(cohort3 = paste0(cohort2, " N=", ct2))

levels <- oo_go %>% gb(cohort3) %>% su(ct = max(n_samples)) %>% ar(desc(ct)) %>% pu(cohort3)
oo_go <- oo_go %>% mu(cohort3 = factor(cohort3, levels = levels))

top_results <- 
oo_go %>% 
 gb(y, cohort3) %>% 
 mu(rk = row_number(desc(abs(z)))) %>% 
 fi(rk <= 5) %>% 
 ug() %>% 
 tm(cohort = cohort3, y, x, covariate, type, 
    est = round(est,2), se = round(se, 2), 
    pval, p_adj_bh = p_adj_fdr, p_adj_by,  samples = ct2, rk) %>% 
 ar(cohort, y)

fwrite(top_results, file = "top_results.csv")

plot <- function( outcome = "Durable Benefit", x_lab = "Log Odds Estimate"){
 oo_go %>% 
  fi(-log10(p_adj_by) < 15) %>% 
  fi(y == outcome) %>% 
  ggplot(aes(x = est, y = -log10(pval), color = type)) + 
  geom_point() +
  facet_wrap(~ cohort3, ncol = 7) +   
  xlim(-2.5,2.5) + 
  xlab(x_lab) +
  ylab("-Log10( P-value )") +
  geom_hline(yintercept = -log10(.001), linetype = "dashed", color = "black") + 
  ggtitle( paste0( outcome )) + 
  go_theme  
}

bor <- plot("Best Response")
dcb <- plot("Durable Benefit")
os <- plot("Overall Survival", x_lab = "Hazard Estimate")
pfs <- plot("PFS Survival", x_lab = "Hazard Estimate")

ds(14, 7)

bor
ggsave("bor_unadjusted.png", width = 14, height = 7)#, dpi = 300)

dcb
ggsave("dcb_unadjusted.png", width = 14, height = 7)

pfs
ggsave("pfs_unadjusted.png", width = 14, height = 7)

os
ggsave("os_unadjusted.png", width = 14, height = 7)

ds(5,5)
oo_go %>% 
 ggplot(aes(x = -log10(pval), y = -log10(p_adj_by))) + geom_point() + 
 go_theme + 
 xlim(0,12) + ylim(0,12) +
 geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "black") + 
 labs(x = "-Log10( P-values )", y = "-Log10( Adjusted P-values )", title = "Adjusted vs Normal P-values")

ds(14, 7)
sig_cohorts <- 
oo_go %>% 
 mu(y = factor(y, levels = c("Best Response","Durable Benefit", "PFS Survival", "Overall Survival"))) %>% 
 gb(y, cohort3) %>% 
 su(sig = mean(p_adj_by < .05)) %>% 
 ar(cohort3) %>% 
 ggplot(aes(x = y, y = sig, fill = y)) + 
 geom_bar(stat = "identity") + 
 facet_wrap(~cohort3, scales = "free_y", ncol = 8) + 
 go_theme + 
 scale_y_continuous(labels = scales::label_percent()) + 
 labs(x = "Outcome Measured", y = "% Significant", title = "Percent Significant by Cohort")

sig_cohorts
ggsave("sig2.png", width = 14, height = 7)

ds(8, 10)

base <- oo_go %>% fi(y == "Best Response", x == "rna_gene_set_tim3" )
base$cohort <- reorder(base$cohort, base$est)

base %>% 
 ggplot(aes(x = est, y = cohort)) +
 geom_point(size = 3, color = "orange") +                # Points for the estimates
 geom_segment(aes(x = est - 1.65*se, xend = est + 1.65*se, y = cohort, yend = cohort), color = "black", size = 1) +  # Confidence intervals
 labs(title = "OS Log Hazard Rate", y = "Cohort", x = "Gene Set Proliferation\nLog Hazard Rate (90% CI)") +
 go_theme + 
 geom_vline(xintercept = 0, col = "red", size = .1)

top_hits_dcb <- 
oo_go %>% 
 fi(y == "Durable Benefit", pval < .05, !grepl("clin", x)) %>% 
 gb(cohort, type, tissue, dna_feature = !grepl("rna_",x)) %>% 
 mu(rk = row_number(pval)) %>% 
 fi(rk <= 5) %>% 
 ug() %>% 
 tm(outcome = y, cohort = cohort, type, tissue, x, dna_feature, est, pval) %>% 
 ar(cohort, desc(dna_feature), pval)

#top_hits_dcb %>% fi(type == "exact")

names(ready)

names(ready$exact$`Anti-PD-1`)
