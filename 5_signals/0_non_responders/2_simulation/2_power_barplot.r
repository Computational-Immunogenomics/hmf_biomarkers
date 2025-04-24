source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

go <- fread("sim_go.csv")

options(repr.plot.width = 10, repr.plot.height = 4)

freq_map <- list("0.05" = "Rare (5%)", "0.1" = "Rare (1%)", "0.3" = "Common (30%)", "0.5" = "Common (50%)" )
sig_map <- list("0" = "None", "-1.3" = "Moderate", "-2.6" = "High", "-3.9" = "Extreme" )

base <- 
go %>% 
 tm(p_lr = Pr...z.., n, beta, b0, p_simple, p_fisher, cov, pp) %>% 
 ga(method, pval, -n, -beta, -b0, -cov, -pp) %>%
 gb(n, beta, b0, cov, method, pp) %>% 
 su(tot = n(), `P-value` = mean(pval < .05), `P-value adjusted` = mean(pval < 7.5014287020248e-05)) %>%
 ug() %>% 
 rw() %>% 
 mu(pp = factor(freq_map[[as.character(pp)]], levels = rev(unlist(unname(freq_map)))), 
    beta = factor(sig_map[[as.character(beta)]], levels = unlist(unname(sig_map))), 
    n = factor(paste0("Sample Size: ", n), levels = paste0("Sample Size: ", n))
   ) %>% 
 ug() %>% 
 mu(  method2 = factor(ifelse(grepl("p_fisher", method), "Fisher's Exact", "Logistic Regression"), 
     levels =  c("Logistic Regression", "Fisher's Exact"))) %>% 
 ga(adjusted, pct_sig, -n, -beta, -b0, -cov, -method, -pp, -tot, -method2) %>% 
 mu(adjusted = factor(adjusted, levels = c( "P-value", "P-value adjusted"))) %>% 
 fi(method2 == "Fisher's Exact")

fills <- c("P-value" = "#7AABD3", "P-value adjusted" = "#e52f28")
alphas <- c("P-value" = .5, "P-value adjusted" = 1)

library(scales) 

p1 <-
base %>% 
 fi(method != "p_simple", !cov, pp %in% c("Rare (5%)", "Common (50%)"), !n %in% c("Sample Size: 50", "Sample Size: 200")) %>% 
 #ggplot(aes(x = as.factor(beta), y = pct_sig, fill = pp, alpha = adjusted)) +
 ggplot(aes(x = as.factor(beta), y = pct_sig, fill = adjusted)) +
 geom_bar(stat = "identity", position = "dodge", color = "black", alpha = .85) + 
 facet_grid(pp ~ n) + 
 scale_fill_manual(values = fills) + 
 #scale_alpha_manual(values = alphas) + 
 go_theme + 
 labs(y = "Statistical Power", x = "Effect Size", title = "Statistical Power Study") + 
 ylim(0,1.1) + 
 scale_y_continuous(labels = percent_format()) +
 geom_text(aes(label = 100*round(pct_sig,2)), vjust = -0.5, position = position_dodge(width = 0.9), size = 3) + 
 coord_cartesian(ylim = c(0, 1.1)) + 
 guides(fill = guide_legend(title = "Method"), alpha = guide_legend(title = "P-value adjustment")) + 
 theme(legend.position = c(.15,.3))

options(repr.plot.width = 7, repr.plot.height = 4)

p1

ggsave( "sim_results.png", plot = p1, width = 7, height = 4)

getwd()
