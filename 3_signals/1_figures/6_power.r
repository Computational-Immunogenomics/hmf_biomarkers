source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(scales)
library(patchwork)

color_map <- c("Non-Response" = "#e52f28", "Response" = "#7AABD3", "20% Baseline Response" = "#e52f28", "40% Baseline Response" = "#7AABD3")
size_map <- list("1%" = .5, "5%" = 1.2, "7%" = 1.5, "10%" = 2, "50%" = 2.5)
alpha_map <- list("P-value signal adjusted" = 1, "Response < 5%" = .4)

extra_theme <- 
theme(axis.text.x = element_text(angle = 0, size = 12, hjust = .5), 
      axis.text.y = element_text(size = 12), 
      plot.title = element_text(size = 16)) 

p_inv <- function(n) (1-(.05)^(1/n))
go <- data.frame( events = seq(5000), p_response_lt = p_inv(seq(5000)))

pts <- data.frame( "x" = c(11,11,11, 59, 149), "y" = c(0,.24, .51,.05, .02))
labels <- data.frame( "x" = c(11,6, 11), "y" = c(-.015,.21, .54))

schema_theme <-  
go_theme + extra_theme + theme(axis.text.x = element_text(hjust = .5, size = 10)) + 
 theme(
  plot.margin = margin(t = 0, r = 60, b = 0, l = 0)  # top, right, bottom, left
)

share <- 
ggplot(go %>% fi(p_response_lt < .55), aes(x = events, y = p_response_lt)) + 
 scale_x_log10(limits = c(4, 800), breaks = c(11, 59, 149, 800), labels = c("11", "59\n108\n890\n5,900", "149\n298\n2,249\n14,900", "Events\nSamples (Prevalence = 50%)\nSamples (Prevalence=7%)\nSamples (Prevalence=1%)\n")) + 
 scale_y_continuous(breaks = c(.02, .05, .24), labels = c("<2%", "<5%", "<24%"), limits = c(-.02,.56)) +
 geom_point() + 
 geom_line() + 
 labs( x = "Number of Events with No Response", 
       y = "Response Probability",
       title = "Non-Response thresholds") + 
 geom_segment(aes(x = 0, xend = 11, y = .24, yend = .24), linetype = "dashed", color = "#e52f28", alpha = .03) +
 geom_segment(aes(x = 11, xend = 11, y = 0, yend = .24), linetype = "dashed", color = "#e52f28", alpha = .3) +
 geom_segment(aes(x = 0, xend = 59, y = .05, yend = .05), linetype = "dashed", color = "grey", alpha = .3) +
 geom_segment(aes(x = 59, xend = 59, y = 0, yend = .05), linetype = "dashed", color = "grey", alpha = .3) + 
 geom_segment(aes(x = 0, xend = 149, y = .02, yend = .02), linetype = "dashed", color = "grey", alpha = .3) +
 geom_segment(aes(x = 149, xend = 149, y = 0, yend = .02), linetype = "dashed", color = "grey", alpha = .3) + 
 schema_theme + 
 geom_point(   data = pts, aes(x = x, y = y), color = "red", size = 3) +
 geom_text(   data = labels, 
           aes(x = x, y = y, label = rev(c("Baseline Response\nNo Event= 51%", "Response\nCI High = 24%", "Estimate = 0% (0/11)"))), 
           color = "#e52f28", size = 3) + 
 geom_text( aes(x = 140, y = .32, label = "Immune Evasion Driver\nICI Therapy Melanoma\n7% Prevalence\n11 non-responder events\n166 Patients", hjust = .5), size = 3, color = "#e52f28") + 
 geom_text( aes(x = 250, y = .18, label = "Upper 95% CI Response Rate\nClopper-Pearson method", hjust = .5), size = 3) + 
 geom_segment(aes(x = 11, y = .51, xend = 59, yend = .32, alpha = .3), color = "#e52f28", alpha = 0.01, size = .1) + 
 geom_segment(aes(x = 11, y = .24, xend = 59, yend = .32), color = "#e52f28", alpha = 0.01, size = .1) + 
 geom_segment(aes(x = 11, y = 0, xend = 59, yend = .32), color = "#e52f28", alpha = 0.01, size = .1) + 
 geom_segment(aes(x = 40, y = .07, xend = 250, yend = .15), color = "black", alpha = 0.01, size = .1) 

go <- 
data.frame( "type" = c("Non-Response", "Non-Response", "Response"),
            event = c(TRUE, FALSE, FALSE), 
            baseline = c(11, 76, 79)) %>% 
 mu("Cohort Size\nNo Response\n<5% Threshold" = baseline * 59/11, 
    "Cohort Size\nNo Response\n<2% Threshold" = baseline * 149/11) %>% 
 ga(threshold, val, -event, -type) %>% 
 mu(threshold = ifelse(threshold == "baseline", "Hartwig Observed", threshold)) %>% 
 mu(Responder = factor(type, levels = rev(c("Non-Response", "Response"))),
    threshold = 
    factor(threshold, levels = c("Hartwig Observed", 
                                 "Cohort Size\nNo Response\n<5% Threshold" , 
                                 "Cohort Size\nNo Response\n<2% Threshold"))) %>%
 mu( `Biomarker Event` = factor(event, levels = c(FALSE, TRUE)))

tots <- go %>% ug() %>% gb(threshold) %>% su(total = round(sum(val))) %>% ug()

colors <- c("Non-Response" = "#e52f28", "Response" = "#7AABD3")

options(repr.plot.height = 6, repr.plot.width = 7.5) 

b <- 
go %>% 
 ggplot(aes(x = threshold, y = val, fill = Responder, alpha = `Biomarker Event`, color = `Biomarker Event`)) + 
 geom_bar(stat = "identity", width = .6) + 
 schema_theme + 
 scale_fill_manual(values = unlist(colors)) + 
 scale_color_manual(values = c("white", "black")) + 
 scale_alpha_manual(values = c(.3, 1)) + 
 labs(x = "Cohort", y = "Patients Observed or Needed", title = "Skin Melanoma ICI Therapy\nHypothetical Required Cohort Sizes") + 
 geom_text(data = tots, aes(x = threshold, label = total, y = total), inherit.aes = FALSE, vjust = -.4) + 
 geom_text( aes(label = round(val)), vjust = 0)

lets_go <- (share | b ) + plot_layout(widths = c(1,1)) 

options(repr.plot.height = 6.5, repr.plot.width = 14) 
share <- 
lets_go + 
plot_annotation(
    title = "Minimum cohort sizes to detect robust non-response biomarker less 5% and 2%",
    subtitle = " ",
    caption = "Assume all patients with non-response marker do not respond, biomarkers is associated with worse odds of response.",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      plot.caption = element_text(size = 10, face = "italic")
    )
  )

ggsave(paste0(FIG_DIR,"minimum_sample_size.png"), plot = share, width = 14, height = 6.5)
