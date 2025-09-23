source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(patchwork)
library(scales)
library(ggcorrplot)
library(survminer)

base <- readRDS(paste0(SHARE_DIR, "3_ready.Rds")) 

fisher_base <- fread(paste0(SHARE_DIR, "fisher_base.csv"))

my_colors <- c("#F04437", "#E81F64", "#903E97", "#65499E", "#4356A5", "#478FCC", "#34A4DD", "#00BCD4", "#009889", "#4BB04F", "#8BC34C", "#CCDA3A", "#FCED3A", "#FFC10E", "#F8991D", "#F1592C", "#7A5649", "#9F9E9E", "#607F8C")

extra_theme <- 
theme(axis.text.x = element_text(angle = 0, size = 12), 
      axis.text.y = element_text(size = 12), 
      plot.title = element_text(size = 16),
      plot.margin = unit(c(1, 1, 1, 0), "cm")) 

alpha_map <- c("TRUE" = 1, "FALSE" = .15)
size_map <- c("TRUE" = 5, "FALSE" = 2, "Maybe" = 2)
shape_map <- c("Worse" = 25, "Better" = 24)
fill_map <- c('Other' = my_colors[18], 'Immunotherapy' = my_colors[1]) 
color_map <- c("FALSE" = "white", "TRUE" = "black")

highlight_interaction_figure <- 
base %>% 
   ar(desc(Treatment)) %>% 
   ggplot( 
    aes(x = prob_response, 
        y = -log10(p_fdr_fisher),
        alpha = focus, 
        color = (p_fdr_surv < .01),
        fill = Treatment, 
        size = focus, 
        shape = Odds)) + 
   geom_point() + 
   scale_alpha_manual(values = alpha_map) + 
   scale_color_manual(values = color_map) +  
   scale_size_manual(values = size_map) + 
   scale_fill_manual(values = fill_map) +  
   scale_shape_manual(values = shape_map) +   
   go_theme + 
   extra_theme + 
   geom_hline(yintercept = 2, alpha = .7, color = "red") + 
   geom_vline(xintercept = .1, alpha = .7, size = .1) + 
   labs(x = "Estimated Probability of Response", 
        y = "Fisher's Test\n-Log10 (p-value FDR Adjusted)", 
        title = "Systematic Analysis") + 
   guides(alpha = "none", 
          size = "none",  
          shape = guide_legend(override.aes = list(size = 4)),
          fill = guide_legend(override.aes = list(shape = 21, size = 4))
         ) + 
   scale_x_continuous( labels = percent_format(accuracy = 1), limits = c(0,1))

options(repr.plot.width = 8, repr.plot.height = 6)
highlight_interaction_figure

cohort_select <- "Pan-Cancer / Immunotherapy"

interaction_base <- 
base %>% 
 drop_na(p_fdr_surv) %>% 
 mu(focus = 
    ifelse(p_fdr_fisher_by < .1 & 
           Odds == "Worse" & 
           p_fdr_surv_by < .1 & 
           cohortGo == cohort_select, 
           TRUE, 
           FALSE), 
    `PFS Significant` = (p_fdr_surv < .1), 
    Treatment = ifelse(focus, mechanism, "Other"))

base %>% fi(grepl("KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", feature), cohortGo == cohort_select)

interaction_base %>% 
 fi(cohortGo == cohort_select, focus, !grepl("GRAFT", feature),
    feature != "rna_geneset_gene_set_cd8_t_effector_lt50", 
    !grepl("GRAFT", feature), 
    !((grepl("lt25", feature) | grepl("lt75", feature)) & grepl("rna", feature)))

highlight_features <- 
interaction_base %>% 
 fi(cohortGo == cohort_select, 
    focus, 
    feature != "rna_geneset_gene_set_cd8_t_effector_lt50", 
    !grepl("GRAFT", feature),
    !((grepl("lt25", feature) | grepl("lt75", feature)) & grepl("rna", feature))) %>% 
 pu(feature)

tmp <-
fisher_base %>% 
 fi(cohortGo == cohort_select) %>% 
 se( any_of(highlight_features) )

replacements <- 
c("_" = " ", 
  "rna geneset " = "RNA", 
  "gene set" = "",
  "HALLMARK" = "", 
  "KEGG" = "",
  "gt0" = "",
  "gt75" = "Very High",
  "gt50" = "High",
  "gt25" = "Mod/High",
  "lt75" = "Low/Mod",
  "lt50" = "Low",
  "lt25" = "Very Low",
  "purity tmbStatus" = "TMB",
  "hotspot KRAS position 25398284" = "KRAS G12D hotspot",
  "purity tmbPerMb lt6" = "TMB per MB < 6",
  "purity tmbPerMb lt8" = "TMB per MB < 8",
  "neo ct" = "RNA Neoantigens",
  "TGF BETA SIGNALING PATHWAY" = "TGFB",
  "purity tmlStatus low" = "TML Low",
  "purity tmbPerMb lt4" = "TMB per MB < 4",
  "RNAAPM" = "RNA APM",
  "low" = "Low",
  "SIGNALING " = "",
  "t cell" = "T-cell",
  "RENIN ANGIOTENSIN SYSTEM " = str_to_title("RENIN ANGIOTENSIN SYSTEM "),
  " AND " = "/",
  "gep" = "GEP",
  "cd8" = "CD8",
  "RNAImm" = "RNA Imm",
  "RNACD" = "RNA CD",
  "CD 8" = "CD8", 
  "BASAL CELL CARCINOMA" = "Basal Cell")

# Loop over names and replace
s1 <- names(tmp)
for (pattern in names(replacements)) {
  s1 <- gsub(pattern, replacements[pattern], s1)
}
names(tmp) <- s1

cor_base <- cor( tmp,use = "pairwise.complete.obs")

options(repr.plot.width = 10, repr.plot.height = 10)
correlations <- 
ggcorrplot(cor_base, 
           method = "circle",   # "circle", "square", "lower", "upper"
           insig = "blank",     # hide insignificant correlations
           tl.cex = 8, 
           hc.order = TRUE) + 
 scale_x_discrete(guide = guide_axis(n.dodge = 5, check.overlap = TRUE)) + 
 theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5), axis.text = element_text(size = 16),  legend.position = "none") +
 ggtitle("Correlation of Highly Significant Features")

options(repr.plot.width = 4, repr.plot.height = 4)
correlations

pfs_base <-
fisher_base %>% 
 fi(cohortGo == cohort_select) %>% 
 tm( sampleId, 
     primaryTumorLocation, 
     daysToPfsEvent, 
     pfsEvent, 
     non_response,  
     tcell = rna_geneset_CD_8_T_effector_lt50,
     tmb = purity_tmlStatus_low,
     hedgehog = rna_geneset_KEGG_BASAL_CELL_CARCINOMA_gt75,
     renin = rna_geneset_KEGG_RENIN_ANGIOTENSIN_SYSTEM_gt50,
     none = as.numeric((tmb+tcell+hedgehog+renin == 0)), 
     tmb_tcell = as.numeric((tmb + tcell == 2) & hedgehog + renin == 0),
     tmb_tcell_renin = as.numeric((tmb + tcell + renin == 3) & (hedgehog == 0)),
     all = as.numeric((tmb+tcell+hedgehog+renin == 4)),
     event = as.factor(tmb+tcell+hedgehog+renin)) %>% 
 mu( location = ifelse(primaryTumorLocation %in% c("Skin", "Lung", "Bladder"), primaryTumorLocation, "Other")) %>%
 drop_na() %>% 
 drop_na(event) %>% 
 gb(event) %>% 
 mu(pct_response = mean(1-non_response)) %>% 
 mu(Events = as.factor(paste0(event, " (", 100*round(pct_response,2), "% response)")))

surv_formula <- expr(Surv(daysToPfsEvent, pfsEvent) ~ Events)
fits <- survfit(eval(surv_formula, envir = environment()), data = pfs_base)
pval <- signif(survdiff(eval(surv_formula, envir = environment()), data = pfs_base)$pvalue,2)

#fit
oo <- 
ggsurvplot(
 fits,
 data = pfs_base,
 #palette = c("#7AABD3", "#e52f28", "#7AABD3", "#e52f28"),
 conf.int = TRUE, 
 risk.table = TRUE, 
 pval.coord = c(700, .95), 
 xlim = c(0, 365*4), 
 break.time.by = 300, 
 ggtheme = theme_minimal(),
 xlab = "Days", 
 ylab = "Survival Probability", 
 title = "Combination Marker\nLow TMB, Low T-cell, High Renin Angiotensin, High Basal Cell" ) 

oo$plot <- 
oo$plot + 
 annotate("text", x = 1000, y = 0.9, label = paste0("Log-rank        \np-value = ", pval), size = 5) + 
 theme(plot.title = element_text(hjust = 0.5), legend.position = "right") + 
 guides(color = guide_legend(nrow = 5, byrow = TRUE)) 

pfs_plot <- oo$plot

los_tres <- (highlight_interaction_figure + theme(legend.position = "none") | 
             correlations + theme(plot.margin = margin(0, 2, 0,0, unit = "cm"), axis.text = element_text(size = 40)) | 
             pfs_plot) + plot_layout(widths = c(1, 1, 1))

los_dos <- (highlight_interaction_figure | correlations ) + plot_layout(widths = c(1, 1))

options(repr.plot.width = 16, repr.plot.height = 6onda )
share <- 
los_tres  + 
plot_annotation(
    title = "Pan-Cancer Immunotherapy Treated (376 Patients)",
    subtitle = "Combine Independent Highly Signficant Markers (TMB, T-cell Infiltration, Hedgehog Signalling, Renin Angiotension)",
    caption = "Full catalogue of output shared with Supplement",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10, face = "italic")
    )
  )

share

ggsave( paste0(FIG_DIR, "combination.png"), plot = share, width = 18, height = 7)
