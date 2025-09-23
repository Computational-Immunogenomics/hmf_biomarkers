source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(patchwork)

go <- 
readRDS(paste0(SHARE_DIR, "3_ready.rds")) %>% 
 fi(selected_example) %>% 
 ar(highlight_fisher, desc(prob_response), e_nr) %>% 
 mu(prevalence = paste0(round(events/total_patients*100), "%")) %>% 
 mu(example_base = paste0(gsub("\n", " (", example), ")"),
    example = paste0(example_base, " (", prevalence, ")"))

examples <- go %>% pu(example)

fwrite( go, paste0(TMP_DIR, "pfs_highlights.csv") )

tmp <- 
go %>%
 fi(example %in% examples) %>% 
 se(example, feature, e_nr, ne_nr, e_r, ne_r, fisher_pval, p_fdr_fisher, p_fdr_surv, surv_pval) %>% 
 ga(event, ct, -feature, -example, -fisher_pval, -p_fdr_fisher, -p_fdr_surv, -surv_pval) %>% 
 mu(events = !m("ne_", event), 
    dcb = !m("_nr", event),
    ct = as.numeric(ct),
    ct_if_dcb = ifelse(dcb, ct, 0), 
    event = factor(event, levels = c( "ne_r", "ne_nr", "e_r", "e_nr")), 
    example = factor(example, levels = examples))

base_non_event <- tmp %>% fi(!events) %>% mu(dcb = factor(dcb, levels = c(TRUE, FALSE)))
base_event <- tmp %>% fi(events) %>% mu(dcb = factor(dcb, levels = c(TRUE, FALSE)))

tot_non_event <- 
base_non_event  %>% 
 gb(example) %>% 
 su(tot_dcb = sum(as.numeric(ct_if_dcb)), 
    tot = sum(as.numeric(ct)), 
    pct_dcb = tot_dcb/tot) %>% ug()

tot_event <- 
base_event  %>% 
 gb(example) %>% 
 su(tot_dcb = sum(as.numeric(ct_if_dcb)), 
    tot = sum(as.numeric(ct)), 
    pct_dcb = tot_dcb/tot, 
    p_fdr = min(p_fdr_fisher), 
    p_fdr_surv = min(p_fdr_surv), 
    p_fisher = min(fisher_pval),
    p_surv = min(surv_pval)) %>% 
 ug()

base_theme <-
theme_minimal() +
 theme(plot.title = element_text(hjust = .5), 
    panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
    axis.title.y = element_blank(), legend.position = "none", 
    axis.text.y = element_text(hjust = .5, size = 11),
    axis.text.x = element_text(hjust = .5, size = 11),
    axis.title.x = element_text(hjust = .5, size = 14)) 

no_y <- theme(axis.text.y = element_blank())
scale_x_event <- scale_x_continuous(expand = expansion(mult = c(0.05, 0.2)), breaks = c(0,25, 50), limits = c(0,63))

alphas <- c("FALSE" = .6, "TRUE" = .8)
colors <- c("FALSE" = "black", "TRUE" = "black")
response <- c("FALSE" = "#e52f28", "TRUE" = "#7AABD3")

settings_stuff <- list(
  scale_fill_manual( values = response),
  scale_alpha_manual( values = alphas),
  scale_color_manual( values = colors), 
  base_theme
)

non_events <- 
base_non_event %>% 
 ggplot( aes( y = example, x = as.numeric(ct), alpha = events, fill = dcb, color = events)) +
 geom_bar(stat = "identity", width = 0.7) + 
 settings_stuff +
 scale_x_reverse(expand = expansion(mult = c(0.25, 0.05))) +
 labs(x = "Number of Patients", title = "No Event") + 
 geom_text(aes(label = ct), position = position_stack(vjust = 0.5),  color = "black") +
 geom_text( data = tot_non_event, aes(label = paste0(round(100*pct_dcb), "%"), x = tot, y = example),  
           color = "black", inherit.aes = FALSE, hjust = 1.3, size = 4) 

events <- 
base_event %>% 
 ggplot( aes(y = example, 
             x = as.numeric(ct), 
             alpha = events, fill = dcb, 
             color = events)) + 
 geom_bar(stat = "identity", width = 0.7) +
 settings_stuff + no_y + scale_x_event +
 labs(x = "Number of Patients", title = "Event") + 
 geom_text(aes(label = ct), position = position_stack(vjust = 0.5),  color = "black") +
 geom_text( data = tot_event, aes(label = paste0(round(100*pct_dcb), "%"), x = tot, y = example),  
           color = "black", inherit.aes = FALSE, hjust = -.4, size = 4) 

options(repr.plot.height = 6, repr.plot.width = 8)

ors <- 
go %>%
 mu(example = factor(example, levels = examples)) %>% 
 ggplot( aes(y = example, x = or)) + 
  geom_point() +                              # points for estimates
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +  # horizontal error bars for CI
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +      # reference line (e.g., OR=1)
  xlab("Odds Ratio Response") +
  ylab("") +
  ggtitle("Odds Ratio 95% CI") + 
  settings_stuff + scale_x_continuous( breaks = c(0,1,2), limits = c(0,2)) + 
  geom_text( data = tot_event, aes(label = paste0("Adj p = ", signif(p_fdr,1), "\nRaw p = ", signif(p_fisher, 1)), 
                                   x = 1.1, y = example), color = "black", inherit.aes = FALSE, hjust = 0, size = 2.5)

hazard <- 
go %>%
 mu(example = factor(example, levels = examples)) %>% 
 ggplot( aes(y = example, x = surv_est)) + 
  geom_point() +                              # points for estimates
  geom_errorbarh(aes(xmin = surv_low, xmax = surv_high), height = 0.2) +  # horizontal error bars for CI
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +      # reference line (e.g., OR=1)
  xlab("Hazard Estimate") +
  ylab("") +
  ggtitle("PFS Hazard 95% CI") + 
  settings_stuff + scale_x_continuous( breaks = c(-2,0,2), limits = c(-3,3)) + 
  geom_text( data = tot_event, 
            aes(label = paste0("Adj p = ", signif(p_fdr_surv,1),"\nRaw p = ", signif(p_surv,1)), 
                x = -2.5, y = example), color = "black", inherit.aes = FALSE, hjust = 0, size = 2.5)

together <- (non_events | events | ors + no_y | hazard + no_y) + plot_layout(widths = c(1, 1, .7,.7)) 

options(repr.plot.height = 7.5, repr.plot.width = 16)
share <- 
together + 
plot_annotation(
    title = "Highlighted Univariate Examples",
    subtitle = "Less than < 5% Estimated Response, only Top 4 Significant for Fisher's Test after FDR Adjustment",
    caption = "Full catalogue of output shared with Supplement",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10, face = "italic")
    )
  )

share

ggsave(paste0(FIG_DIR, "highlights.png"), share, width = 16, height = 7.5)
