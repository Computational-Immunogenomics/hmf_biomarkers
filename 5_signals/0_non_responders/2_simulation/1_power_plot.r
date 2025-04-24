source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(scales)

go <- fread("sim_go_new.csv")

head(go)

#go %>% gb(n, p_base, p_event) %>% su(ct = n()) %>% ar(desc(ct))

plts_base <- 
go %>% 
 #fi(p_event == 0) %>% 
 mu(`  Feature\nPrevalence` = as.factor(prevalence)) %>% 
 gb(`  Feature\nPrevalence`, n, p_base, p_event) %>% 
 su( raw = mean(p_fisher < .05), adjusted = mean(p_fisher < .001), never_response = mean(ci.high < .02, na.rm = TRUE))

mapper <- c("0" = "Probability Response = 0% (Never Response)",
  "0.01" = "Probability Response = 1%",
  "0.1" = "Probability Response = 10%",
  "0.2" = "Probability Response = 20%",
  "0.4" = "Probability Response = 40% (No Signal)")

mapper_threshold <- 
c("raw" = "P-value raw", 
  "adjusted" = "P-value adjusted", 
  "never_response" = "Response < 2%")

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

options(repr.plot.width = 6, repr.plot.height = 4)

plts_ready %>% 
 fi(p_event == 0) %>% 
 ggplot( aes(x = n, y = val, alpha = threshold, color = `  Feature\nPrevalence`)) + 
 geom_point(size = 3) +
 geom_line(aes(group = interaction(threshold,`  Feature\nPrevalence`)), linewidth = 1.2) + 
 go_theme + 
 scale_x_continuous(trans = "log10", breaks = c(20, 40, 80, 100, 200, 500, 1000, 2000, 5000)) + 
 labs(y = "Statistical Power",
      x = "Sample Size", 
      title = "Statistical Power to detect never response signals") +  
 scale_y_continuous(labels = label_percent()) #+ scale_color_brewer(palette = "Set4")

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

plts_ready %>% 
 rw() %>% mu(event = factor(mapper[as.character(p_event)], levels = rev(unname(mapper)))) %>% ug() %>%
 ggplot( aes(x = expected_events, y = val, alpha = gp, color = `  Feature\nPrevalence`)) + 
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

plts_ready %>% 
 rw() %>% mu(event = factor(mapper[as.character(p_event)], levels = rev(unname(mapper)))) %>% ug() %>%
 ggplot( aes(x = expected_non_events, y = val, alpha = threshold, color = `  Feature\nPrevalence`)) + 
 geom_point(size = 3) +
 geom_line(aes(group = interaction(threshold,`  Feature\nPrevalence`)), linewidth = 1.2) + 
 facet_wrap(~event, ncol = 5) + 
 go_theme + 
 scale_x_continuous(trans = "log10", breaks = c(1, 5, 10, 20, 40, 100, 500, 1000), limits = c(1,500)) + 
 labs(y = "Statistical Power",
      x = "Expected Number of Non-Events = (Sample Size * (1-Feature Prevalence))", 
      title = "Statistical Power to detect response signals") +  
 scale_y_continuous(labels = label_percent()) + 
 geom_hline(yintercept = .05) + 
 geom_text(data = annotate, aes( label = name), alpha = 1, color = "black", size = 3)

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
