source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "fisher.r"))

library(survminer)
library(rlang)
library(cowplot)
library(patchwork)
library(ggpubr)

go <- fread(paste0(SHARE_DIR, "fisher_base.csv"))

examples <- fread(paste0(TMP_DIR, "pfs_highlights.csv")) %>% se(cohortGo, feature, example_base)

ready_data <- function( i, j ) {
 go %>% 
  fi(cohortGo == i) %>% 
  se(contains("Pfs"), !!sym(j)) %>% 
  rename( val = !!j) %>% 
  mu(cohortGo = !!i, feature = !!j)
}

base <- data.frame()
for( i in seq(nrow(examples))){
  tmp <- ready_data(examples[i,]$cohortGo, examples[i,]$feature) %>% drop_na()
  base <- rbind(base, tmp)
}

ready <- base %>% lj(examples, by = c("cohortGo", "feature"))

data_prep <- function(example) {
    data <- ready %>% fi( example_base == !! example) 
    surv_formula <- expr(Surv(daysToPfsEvent, pfsEvent) ~ val)
    fits <- survfit(eval(surv_formula, envir = environment()), data = data)
    pval <- signif(survdiff(eval(surv_formula, envir = environment()), data = data)$pvalue,2)
    list("data" = data, "fits" = fits, "surv_formula" = surv_formula, "pval" = pval)
}

figurer <- function(example, data, fits, pval, pad = 30 ) {
        
    max_time <- max(data$daysToPfsEvent, na.rm = TRUE)
    xmax <- min(max_time, 1100)
    
    oo <- 
    ggsurvplot(
     fits,
     data = data,
     palette = c("#7AABD3", "#e52f28"),
     conf.int = TRUE, 
     risk.table = TRUE, 
     pval.coord = c(700, .95), 
     xlim = c(0, xmax), 
     break.time.by = 300, 
     ggtheme = theme_minimal(),
     xlab = "Days", 
     ylab = "Survival Probability", 
     title = example) 

    oo$plot <- 
    oo$plot + 
     annotate("text", x = 500, y = 0.9, label = paste0("Log-rank        \np-value = ", pval), size = 5) + 
     theme(plot.title = element_text(hjust = 0.5)) + 
     guides(color = guide_legend(nrow = 1, byrow = TRUE)) 

    as_ggplot(gridExtra::arrangeGrob(oo$plot, oo$table, layout_matrix = matrix(c(1,1,1,1,1,1,1,2,2)))) + 
    theme(plot.margin = margin(t = pad, r = pad, b = pad, l = pad))
}

survival_figures <- list()
for( i in rev(unique(ready$example_base))){
    data_ready <- data_prep(i)
    data <- data_ready$data 
    fits <- data_ready$fits
    pval <- data_ready$pval
    surv_formula <- data_ready$surv_formula
    survival_figures[[i]] <- figurer( i, data = data, fits = fits, pval = pval)
}

options(repr.plot.height = 24, repr.plot.width = 24) 
wrap_plots(survival_figures, ncol = 4)

ggsave(paste0(FIG_DIR, "highlights_pfs.png"), wrap_plots(survival_figures, ncol = 4), width = 24, height = 24)
