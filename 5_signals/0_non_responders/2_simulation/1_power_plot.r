source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(scales)

go <- fread("sim_go_new.csv")

go_binom_test <- function( n, x, p = .02) {
   if(is.na(x)){1}
   else if (n == 0) {1}
   else{ binom.test(x, n, p, alternative = "less")$p.value } 
}

binom.test(0, 1, p = 0.05, alternative = "less")$p.value

tail(go)

plts_base <- 
go %>% 
 #fi(p_event == 0) %>% 
 mu( prevalence = as.factor(prevalence), `  Feature\nPrevalence` = prevalence) %>% 
 gb(`  Feature\nPrevalence`, n, p_base, p_event) %>% 
 su( signal_raw = mean(p_fisher < .05), 
     signal_adjusted = mean(p_fisher < .004),
     never_response_lt_10 = mean(ci.high < .1, na.rm = TRUE),
     never_response_lt_05 = mean(ci.high < .05, na.rm = TRUE),
     never_response_lt_02 = mean(ci.high < .02, na.rm = TRUE))

mapper <- c(
  "0" = "Probability Response = 0% (Never Response)",
  "0.01" = "Probability Response = 1%",
  "0.1" = "Probability Response = 10%",
  "0.2" = "Probability Response = 20%",
  "0.4" = "Probability Response = 40% (No Signal)")

mapper_threshold <- 
c("signal_raw" = "P-value signal raw", 
  "signal_adjusted" = "P-value signal adjusted", 
  "never_response_lt_10" = "Response < 10%",
  "never_response_lt_05" = "Response < 5%",
  "never_response_lt_02" = "Response < 2%")

plts_ready <- 
plts_base %>% 
 ga(threshold, val, -`  Feature\nPrevalence`, -n, -p_base, -p_event) %>% 
 mu(expected_non_events = n * (1-as.numeric(as.character(`  Feature\nPrevalence`))), 
    expected_events = n * as.numeric(as.character(`  Feature\nPrevalence`)), 
    expected_events_non_response = n * as.numeric(as.character(`  Feature\nPrevalence`)) * p_event, 
    expected_events_response = n * as.numeric(as.character(`  Feature\nPrevalence`)) * p_base) %>% 
 rw() %>% 
 mu(event = factor(mapper[as.character(p_event)], levels = rev(unname(mapper))), 
    gp = factor(mapper_threshold[[threshold]], levels = rev(unname(mapper_threshold)))) %>% ug()

options(repr.plot.width = 8, repr.plot.height = 4)

p1 <- 
plts_ready %>% 
 fi(p_event == 0, !grepl("never", threshold)) %>% 
 ggplot( aes(x = n, y = val, alpha = gp, color = `  Feature\nPrevalence`)) + 
 geom_point(size = 3) +
 geom_line(aes(group = interaction(gp,`  Feature\nPrevalence`)), linewidth = 1.2) + 
 go_theme + 
 scale_x_continuous(trans = "log10", breaks = c(20, 30, 40, 50, 60, 80, 100, 200, 500, 1000, 2000, 5000)) + 
 labs(y = "Statistical Power",
      x = "Sample Size", 
      title = "Statistical Power to detect Never Response signals") +  
 scale_y_continuous(labels = label_percent()) 

s1 <- 
plts_ready %>% 
 fi(val > .8) %>% 
 gb(gp, `  Feature\nPrevalence`) %>% 
 su( min_samples = min(n), .groups = "drop") %>% 
 ug() 

power_summary <- 
s1 %>% 
 complete(gp, `  Feature\nPrevalence`, fill = list(min_samples = 5000)) %>% 
 mu(label = ifelse(min_samples == 5000, "5000+", as.character(min_samples)))

p2 <- 
power_summary %>% 
 ggplot( aes(x = `  Feature\nPrevalence`, y = min_samples, fill = gp)) + 
 geom_bar(stat = "identity", position = "dodge", color = "black") + 
 geom_text(aes(label = label), position = position_dodge(width = 0.9),vjust = -0.5, size = 4) + 
 go_theme + 
 ylim(0, 5200) + 
 labs( x = "Biomarker Prevalence", y = "# Patients for 80% Power", title = "Biomarkers for Never Response\n(40% Baseline Response Rate)")

options(repr.plot.width = 10, repr.plot.height = 4)
p2

options(repr.plot.width = 16, repr.plot.height = 4)

annotate <- 
rbind(
    plts_ready %>% fi(expected_events == 30, event == "Probability Response = 10%", n == 100) %>% 
    mu(name = "TMB Low\n(NSCLC Anti-PD1)"), 
    plts_ready %>% fi(expected_events == 10, event == "Probability Response = 0% (Never Response)", n == 100) %>% 
    mu(name = "B2M Loss\n(Melanoma Anti-PD1)") , 
    plts_ready %>% fi(event == "Probability Response = 20%", n == 500, expected_events == 150) %>% 
    mu(name = "TMB Low\n(Anti-PD1 Pan-Cancer)")  
)    

p3 <- 
plts_ready %>% 
 fi(`  Feature\nPrevalence` == .3) %>% 
 rw() %>% mu(event = factor(mapper[as.character(p_event)], levels = rev(unname(mapper)))) %>% ug() %>%
# ggplot( aes(x = expected_events, y = val, alpha = gp, color = `  Feature\nPrevalence`)) + 
 ggplot( aes(x = expected_events, y = val, color = gp)) + 
 geom_point(size = 3) +
 geom_line(aes(group = interaction(threshold,`  Feature\nPrevalence`)), linewidth = 1.2) + 
 facet_wrap(~event, ncol = 5) + 
 go_theme + 
 scale_x_continuous(trans = "log10", breaks = c(1, 5, 10, 20, 40, 100, 500, 1000), limits = c(1,500)) + 
 labs(y = "Statistical Power",
      x = "Expected Number of Events = (Sample Size * Feature Prevalence)", 
      title = "Statistical Power to detect response signals") +  
 scale_y_continuous(labels = label_percent()) + 
 geom_hline(yintercept = .05) + 
 geom_text(data = annotate, aes( label = name), alpha = 1, color = "black", size = 3)

ggsave("p3.png", plot = p3, width = 16, height = 4)

plts_ready %>% 
 gb(gp, expected_events, event) %>% 
 su(val = mean(val)) %>% 
 ggplot( aes(x = expected_events, y = val, color = gp)) + 
 geom_point(size = 3) +
 geom_line(aes(group = gp), linewidth = 1.2) + 
 facet_wrap(~event, ncol = 5) + 
 go_theme + 
 scale_x_continuous(trans = "log10", breaks = c(1, 5, 10, 20, 40, 100, 500, 1000), limits = c(1,1000)) + 
 labs(y = "Statistical Power",
      x = "Expected Number of Events = (Sample Size * Feature Prevalence)", 
      title = "Statistical Power to detect response signals") +  
 scale_y_continuous(labels = label_percent()) + 
 geom_hline(yintercept = .05) + 
 geom_text(data = annotate, aes( label = name), alpha = 1, color = "black", size = 3) + 
   geom_rect(data = highlight_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightgreen", alpha = 0.3,
            inherit.aes = FALSE) 

plts_ready <- 
plts_base %>% 
 ga(threshold, val, -`  Feature\nPrevalence`, -n,  -p_base, -p_event) %>% 
 mu(expected_events = n * as.numeric(as.character(`  Feature\nPrevalence`)), 
    expected_events_non_response = n * as.numeric(as.character(`  Feature\nPrevalence`)) * p_event, 
    expected_events_response = n * as.numeric(as.character(`  Feature\nPrevalence`)) * p_base) %>% 
 rw() %>% mu(event = factor(mapper[as.character(p_event)], levels = rev(unname(mapper)))) %>% ug()

#plts_ready
