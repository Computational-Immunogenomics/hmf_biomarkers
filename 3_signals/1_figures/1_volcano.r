source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(patchwork)
library(scales)
library(ggrepel)
library(ggcorrplot)
library(corrplot)
library(survminer)

extra_theme <- 
theme(axis.text.x = element_text(angle = 0, size = 12), 
      axis.text.y = element_text(size = 12), 
      plot.title = element_text(size = 16),
      plot.margin = unit(c(1, 1, 1, 0), "cm")) 

base <- readRDS(paste0(SHARE_DIR, "3_ready.rds"))

cohorts <- fread(paste0(SHARE_DIR, "top_mechanisms.csv")) 

options(repr.plot.height = 1, repr.plot.width = 6) 
my_colors <- c("#F04437", "#E81F64", "#903E97", "#65499E", "#4356A5", "#478FCC", "#34A4DD", "#00BCD4", "#009889", "#4BB04F", "#8BC34C", "#CCDA3A", "#FCED3A", "#FFC10E", "#F8991D", "#F1592C", "#7A5649", "#9F9E9E", "#607F8C")
df <- data.frame(color = my_colors,x = seq_along(my_colors))

alpha_map <- list("Both Signficant" = 1, "Fisher Signficant" = .8, "PFS Signficant" = .65, "Both Significant / Unadjusted" = .5, "Rest" = .25)
size_map <- list("Both Signficant" = 5.5, "Fisher Signficant" = 4.5, "PFS Signficant" = 3.5, "Both Significant / Unadjusted" = 2.5, "Rest" = 1.5)
size_map2 <- list("Both Signficant" = 11, "Fisher Signficant" = 9, "PFS Signficant" = 7, "Both Significant / Unadjusted" = 5, "Rest" = 3)
shape_map <- list("Worse" = 25, "Better" = 24)

fill_map <- 
list(
' ' = my_colors[18],
'Other' = my_colors[18],
'none' = my_colors[18],
'significant' = my_colors[18],
'Anti-PD-1' = my_colors[1], 
'Anti-CTLA-4 / Anti-PD-1' = my_colors[1],    
'Immunotherapy' = my_colors[1], 
'Multikinase inhibitor' = my_colors[10],
'Anti-AR' = my_colors[3],     
'Anti-EGFR' = my_colors[4], 
'Aromatase inhibitor' = my_colors[5],
'Hormonal therapy' = my_colors[5],    
'CAPEOX + Bevacizumab'  = my_colors[6],
'Anti-VEGF / Platinum/Pyrimidine (ant)agonist'  = my_colors[7],   
'Alkaloid / Platinum' = my_colors[8],
'Alkylating / Anthracycline' = my_colors[8],    
'Platinum / Pyrimidine (ant)agonist'= my_colors[9],
'Platinum / Taxane'= my_colors[10],
'Pyrimidine (ant)agonist'= my_colors[11],
'Taxane'= my_colors[12],
'Chemotherapy' = my_colors[8], 
'Targeted therapy' = my_colors[10]) 

color_map <- list()
for( i in names(fill_map)){
 if( i %in% c("other", "none", "significant")){
    color_map[[i]] <- my_colors[18]
 } else {
    color_map[[i]] <- "black"
 }
}

main_figure <- 
 base %>% 
   ar(desc(Treatment)) %>% 
   ggplot( 
    aes(x = prob_response, 
        y = -log10(p_fdr_fisher),
        alpha = gp, 
        fill = Treatment, 
        size = gp, 
        shape = Odds
       )) + 
   geom_point(color = "grey") + 
   scale_alpha_manual(values = unlist(alpha_map)) + 
   scale_size_manual(values = unlist(size_map)) + 
   scale_fill_manual(values = unlist(fill_map)) +  
   scale_shape_manual(values = unlist(shape_map)) +   
   go_theme + 
   extra_theme + 
   geom_hline(yintercept = 1, alpha = .7, color = "red") + 
   geom_vline(xintercept = .05, alpha = .7, size = .1) + 
   geom_vline(xintercept = .1, alpha = .2, size = .1) + 
   labs(x = "Estimated Probability of Response", y = "Fisher's Test\n-Log10 (FDR Adjusted p-value)", title = "Systematic Analysis Results") + 
   guides(alpha = "none", 
          color = "none",
          shape = guide_legend(override.aes = list(size = 4)),
          fill = guide_legend(override.aes = list(shape = 21, size = 4))
         ) + 
   scale_x_continuous( labels = percent_format(accuracy = 1), limits = c(0,1))

options(repr.plot.width = 8, repr.plot.height = 6)
main_figure

remove <- c("Pan-Cancer / Pazopanib\nRNA Arachidonic Acid Metabolism Very Low",
            "Lung NSCLC / Chemotherapy\nDrivers Pathway DDR")

zoom_df <- base %>% fi(selected_example) %>% ar(gp) %>% fi(!example %in% remove)

highlight <- 
main_figure + 
 scale_x_continuous( labels = percent_format(accuracy = 1), breaks = c(0,.05), limits = c(-.04,.051))  + 
 ylim(.35,1.5) + 
 scale_size_manual(values = unlist(size_map2)) + 
 theme(axis.title.y = element_blank()) + 
 labs(x = "Estimated Probability of Response", 
      y = "Fisher's Test of Odds Ratio\n-Log10 (p-value)", 
      title = "Highlighted Univariate Results", size = NULL, alpha = NULL, color = NULL) +
 geom_point() + 
 geom_text_repel(data = zoom_df, 
                 aes(label = example), size = 2.5,  nudge_y = .1,
                    force = 1.5,               # increase repulsion force (default = 1)
                    force_pull = 0.1,        # reduce attraction toward anchor point (default = 0.1)
                    box.padding = 0.5, 
                 max.overlaps = Inf) 

go <- (main_figure + theme(legend.position = "none") | highlight ) + plot_layout(widths = c(6, 6))

options(repr.plot.width = 14, repr.plot.height = 7)
share <- go + 
plot_annotation(
    title = "Non-Response - Systematic Testing Results",
    subtitle = "2,561 biomarkers tested for Non-Response association (Fisher's Exact, Cox-PH) across across 55 cohorts",
    #caption = "Highlighted Examples Selected Based on Signficance and Manual Selection",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10, face = "italic")
    )
  )
share

ggsave( paste0(FIG_DIR, "volcano_main.png"), plot = share, width = 14, height = 7)
