source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(patchwork)

top_mechanisms <- fread(paste0(SHARE_DIR, "top_mechanisms.csv"))

top_pan <- top_mechanisms %>% fi(grepl("Pan-Cancer", cohortGo)) 
top_tissue <- top_mechanisms %>% fi(grepl("Pan-Cancer", cohortGo)) 

p_base <- 
top_mechanisms %>% 
 mu(cohortGoo = ifelse(cohortGo == "Pan-Cancer", "All treatments", cohortGo)) %>% 
 mu(pan = grepl("Pan-Cancer", cohortGo) | cohortGo == "All treatments", 
    cohortGoo = gsub("Pan-Cancer ## ", "", cohortGoo), 
    cohortGoo = gsub("\\(ant\\)agonist", "", cohortGoo), 
    cohortGoo = gsub(" ## ", " / ", cohortGoo)) %>% 
 mu(cohortGooo = ifelse(!pan, sub(" / ", "\n", cohortGoo), cohortGoo),
    cohortGoo = gsub("Pan-Cancer / ", "", cohortGoo))

fwrite( p_base , paste0(SHARE_DIR, "cohort_names_map.csv"))

base <- 
p_base %>% 
 ga( response, val, -cohortGo, -group, -ct, -pan, -cohortGoo, -cohortGooo) %>% 
 mu(response = factor(ifelse(response == "no_dcb", "No DCB", "DCB"), levels = rev(c("DCB", "No DCB"))),
    cohortGo = factor(cohortGo, levels = p_base$cohortGo[order(p_base$ct)]), 
    cohortGoo = factor(cohortGoo, levels = p_base$cohortGoo[order(p_base$ct)]),
    group = factor(group, levels = c("treatment", "mechanism", "type"))) %>% 
 group_by(cohortGo, cohortGoo) %>%
 arrange(cohortGo, rev(response)) %>%
 mutate(pos = cumsum(val) - 0.5 * val)

colors <- c("No DCB" = "#e52f28", "DCB" = "#7AABD3")

extra_theme <- 
theme(legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12), 
    axis.text.x = element_text(size = 12, angle = 0), 
    plot.title = element_text(size = 16) ) 

plter <- function(d_base, title = "X", text_size = 4){
 df_totals <- d_base %>% gb(cohortGoo) %>% su(total = sum(val)) 
 ggplot(d_base, aes(x = val, y = cohortGoo, fill = response)) +
  geom_bar(stat = "identity", color = "black", alpha = .8) +
  labs(title = title) +
  geom_text(aes(label = val, x= pos), color = "black", position = position_dodge(), size = text_size) + 
  scale_fill_manual( values = colors) +  
  geom_text(data = df_totals, aes(y = cohortGoo, x = total, label = total), inherit.aes = FALSE, hjust = -.2, size = 4) + 
  scale_x_continuous(expand = expansion(c(0, 0.12))) + 
  go_theme + extra_theme
}

p1 <- plter(base %>% fi(!pan, group == "treatment"), "Cohort Specific: Drug Treatment")
p2 <- plter(base %>% fi(!pan, group == "mechanism"), "Cohort Specific: Drug Mechanism") 
p3 <- plter(base %>% fi(!pan, group == "type"), "Cohort Specific: Drug Type") + labs(title = "Cohort Specific: Drug Type", x = "Number of Patients") + go_theme + extra_theme + theme(legend.position = "none",axis.title.y = element_blank()) 
p4 <- plter(base %>% fi(pan, group == "treatment"), "Pan-Cancer: Drug Treatment") + theme(legend.position = c(0.7, 0.25) ) + labs(fill = NULL)
p5 <- plter(base %>% fi(pan, group == "mechanism"), "Pan-Cancer: Drug Mechanism") 
p6 <- plter(base %>% fi(pan, group == "type"), "Pan-Cancer: Drug Type")  + labs(title = "Pan-Cancer: Drug Type", x = "Number of Patients") + go_theme + extra_theme + theme(legend.position = "none",axis.title.y = element_blank()) + scale_x_continuous(expand = expansion(c(0, 0.2))) 

options(repr.plot.width = 14, repr.plot.height = 8)

cohort <- (p1 / p2 / p3 ) + plot_layout(heights = c(6, 6, 7)) 
pan_cancer <- (p4 / p5 / p6 ) + plot_layout(heights = c(6, 6, 2.5)) 
cohorts_full = (cohort | pan_cancer ) + plot_layout(widths = c(1,1)) 

treatment <- (p1 / p4 ) + plot_layout(heights = c(3,5)) 
mechanism <- (p2 / p5) + plot_layout(heights = c(3,3)) 

cohorts_reduced = (treatment | mechanism ) + plot_layout(widths = c(1,1)) 

options(repr.plot.width = 14, repr.plot.height = 11)
#cohorts_reduced
share <- 
cohorts_full + 
plot_annotation(
    title = "Cohort Sizes of Patient Responders vs Non-Responders",
    subtitle = "Responder status defined with RECIST (CR, PR, SD for 6+ months)",
    caption = "Cohorts required to have 30 total patients and 15 responder and 15 non-responders",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10, face = "italic")
    )
  )

share

ggsave( paste0(SHARE_DIR, "0_cohorts_full.png"), plot = share, width = 14, height = 11)
