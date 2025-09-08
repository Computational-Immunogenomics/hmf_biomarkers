source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "fisher.r"))

library(patchwork)

extra_theme <- 
theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) + 
  theme(axis.text.y = element_text(size = 12), plot.title = element_text(size = 13)) 

ready <- readRDS(paste0(SHARE_DIR, "3_ready.rds"))

base <- 
ready %>% 
 mu(cohortGo = ifelse(cohortGo != "Pan-Cancer", gsub(" \\(ant\\)agonist", "", cohortGo), "All treatments")) %>% 
 gb(cohortGo, pan, group) %>% 
 su( patients = max(total_patients, na.rm = TRUE),
     fisher = sum(sig_fisher & low_responder, na.rm = TRUE), 
     pfs = sum(sig_pfs & low_responder, na.rm = TRUE), 
     `Both P-Adj < .1` = sum(sig_fisher + sig_pfs + low_responder == 3)) %>% 
 mu( `Fisher Only P-Adj < .1` = fisher - `Both P-Adj < .1`,
     `Cox PH Only P-Adj < .1` = pfs - `Both P-Adj < .1`) %>% 
 se( -fisher, -pfs) %>% 
 ga(gp, val, -cohortGo, -group, -pan, -patients) %>% 
 ug() %>%
 mu( gp = factor(gp, levels = rev(c("Cox PH Only P-Adj < .1","Fisher Only P-Adj < .1", "Both P-Adj < .1"))))

base %>% gb(gp) %>% su(sum(val))

plter <- function(plt_base, title = "X", text_size = 4){
 dd_tot <- plt_base %>% gb(cohortGo, patients) %>% su(tot = sum(val))
 plt_base %>%
  ggplot(aes(y = fct_reorder(cohortGo, patients), x = val, alpha = gp)) + 
  geom_bar(stat = "identity", color = "black", fill = "#e52f28") +   
  go_theme + 
  extra_theme + 
  labs(title = title, alpha = NULL) + 
  scale_x_continuous(expand = expansion(c(0, .05)), limits = c(0, 7), breaks = seq(7)) +  
  geom_text(data = dd_tot, aes(y = cohortGo, x = tot, label = tot), inherit.aes = FALSE, hjust = -.2, size = 4) + 
  scale_alpha_manual(values = c("Both P-Adj < .1" = 1, "Fisher Only P-Adj < .1" = .6,  "Cox PH Only P-Adj < .1" = .2)) 
}

p1 <- plter(base %>% fi(!pan, group == "treatment"), "Cohort Specific: Drug Treatment") 
p2 <- plter(base %>% fi(pan, group == "treatment"), "Pan-Cancer: Drug Treatment") + theme(legend.position = c(0.5, 0.3), legend.text = element_text(size = 13) ) + labs(fill = NULL)
p3 <- plter(base %>% fi(!pan, group == "mechanism"), "Cohort Specific: Drug Mechanism")
p4 <- plter(base %>% fi(pan, group == "mechanism"), "Pan-Cancer: Drug Mechanism")
p5 <- plter(base %>% fi(!pan, group == "type"), "Cohort Specific: Treatment Type") 
p6 <- plter(base %>% fi(pan, group == "type"), "Pan-Cancer: Treatment Type") 

options(repr.plot.width = 14, repr.plot.height = 8)

cohort <- (p1 / p3 / p5 ) + plot_layout(heights = c(6, 6, 7)) 
pan_cancer <- (p2 / p4 / p6 ) + plot_layout(heights = c(6, 6, 2.5)) 
cohorts_full = (cohort | pan_cancer ) + plot_layout(widths = c(1,1)) 

options(repr.plot.width = 14, repr.plot.height = 11)
share <- 
cohorts_full + 
plot_annotation(
    title = "Number of Significant Non-Response Markers after FDR Adjustment",
    subtitle = "Less than < 5% Estimated Response (2,561 markers tested across cohorts)",
    caption = "FDR applied separately for PFS and Fisher's test across whole analysis",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10, face = "italic")
    )
  )

share 

ggsave( paste0(FIG_DIR, "cohort_summaries.png"), plot = share, width = 14, height = 11)
