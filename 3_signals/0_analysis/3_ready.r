source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))
source(paste0(HELP_DIR, "names.r"))

base <- fread(paste0(SHARE_DIR, "2_run_marginal_output.csv")) 

cohorts <- fread(paste0(SHARE_DIR, "top_mechanisms.csv"))

treatment_mechanism_map <- 
fread(paste0(SHARE_DIR, "treatment_mechanism_map.csv")) %>% 
 mu(derived_treatmentName = gsub( "##","/", derived_treatmentName),
    derived_treatmentMechanism = gsub( "##","/", derived_treatmentMechanism))

low_responder_threshold <- .05
pval_threshold <- .1

orderer <- c('Anti-PD-1','Immunotherapy','Chemotherapy','Anti-AR',' ')

lets_go <- 
base %>% 
 fi(total_patients >= 30, non_responders >= 15, responders >= 15, events >= 5) %>% 
 mu(p_fdr_fisher = p.adjust(fisher_pval, method = "fdr"), 
    p_fdr_fisher_by = p.adjust(fisher_pval, method = "BY"),
    p_fdr_surv = p.adjust(surv_pval, method = "fdr"),
    p_fdr_surv_by = p.adjust(surv_pval, method = "BY"),
    or = ifelse(or == "Inf", exp(5), or), 
    surv_high = surv_est + 1.96*surv_se, surv_low = surv_est - 1.96*surv_se,
    prob_response = e_r/events,
    low_responder = (prob_response <= low_responder_threshold),
    sig_fisher = (p_fdr_fisher <= pval_threshold), sig_pfs = (p_fdr_surv <= pval_threshold),
    light_highlight = ((low_responder) & (fisher_pval < .02) & (surv_pval < .05)),
    highlight_fisher = ((low_responder) & (sig_fisher) & (surv_pval < .05)), 
    highlight_pfs = ((low_responder) & (sig_pfs) & (fisher_pval < .05)),
    dark_highlight = ((highlight_fisher ) & (highlight_pfs)),
    Odds = ifelse(direction == "Response", "Better", "Worse"),
    pan = grepl("Pan-Cancer", cohortGo)) %>% 
  fi(surv_se <= 10) %>% 
  rw() %>% mu( derived_treatmentName = str_split_fixed( cohortGo	, " / ", n = 2)[2]) %>% ug() %>% 
  lj( cohorts %>% se(cohortGo, group), by = "cohortGo" ) %>% 
  lj( treatment_mechanism_map , by = "derived_treatmentName") %>% 
  mu( derived_treatmentMechanism = ifelse(is.na(derived_treatmentMechanism), derived_treatmentName, derived_treatmentMechanism),
      treatment = derived_treatmentName, 
      mechanism = gsub(" ## ", "/", derived_treatmentMechanism),
      Treatment = factor(ifelse(highlight_fisher, mechanism, " "), levels = orderer)  ) %>% 
  se(-direction) %>% 
  mu(
  gp = factor(
    case_when(
     dark_highlight ~ "Both Signficant",
     highlight_fisher ~ "Fisher Signficant",  
     highlight_pfs ~ "PFS Signficant",   
     light_highlight ~ "Both Significant / Unadjusted",
     TRUE ~ "Rest"), 
       levels = c("Both Signficant","Fisher Signficant","PFS Signficant", "Both Significant / Unadjusted","Rest")))

dim(lets_go)

lets_go %>% ar(fisher_pval) %>% fi(low_responder)

tmp <- 
lets_go %>% 
 fi(highlight_fisher | highlight_pfs | (grepl("driver", feature) & light_highlight)) %>% 
 fi( !((cohortGo == "Skin Melanoma / Anti-PD-1") & (feature %in% c("drivers_pathway_IMMUNE_EVASION_gt0", "drivers_immune_evasion_gt0")))) %>%
 gb(cohortGo) %>%
 mu(rk = row_number(desc(e_nr))) %>%
 tm(cohortGo, feature, fisher_pval, surv_pval, e_nr, e_r, select_example = TRUE, rk) %>% 
 ar(cohortGo, fisher_pval)

remove_examples <- c("signature_log_ID2_gt75", "lilac_hla_cn_B44_gt0", "rna_geneset_KEGG_CYSTEINE_AND_METHIONINE_METABOLISM_gt75")

examples_base <- tmp %>% fi(rk <= 1, !feature %in% remove_examples) %>% tm(cohortGo, feature, selected_example = TRUE, rk)

examples <- 
examples_base %>%
 bind_rows(
 lets_go %>% 
  fi(cohortGo == "Prostate / Anti-AR", feature == "rna_geneset_KEGG_TGF_BETA_SIGNALING_PATHWAY_gt50") %>% 
   tm(cohortGo, feature, selected_example = TRUE) 
 ) %>% 
 bind_rows(
 lets_go %>% 
  fi(cohortGo == "Colorectum / Chemotherapy", feature == "hotspot_KRAS_G12D") %>% 
   tm(cohortGo, feature, selected_example = TRUE) 
 )

lets_go_with_examples <- lets_go %>% lj(examples, by = c("cohortGo", "feature"))

highlights <- examples %>% pu(feature)

s1 <- highlights
for (i in names(update_names)) {
  s1 <- gsub(i, update_names[i], s1)
}

s1 <- str_to_title(s1)
for (i in names(update_names_again)) {
  s1 <- trimws(gsub(i, update_names_again[i], s1))
}

highlights_go <- setNames(s1, highlights)

next_step <- 
lets_go_with_examples %>% 
 rw() %>% 
 mu(clean_name = ifelse(feature %in% names(highlights_go), highlights_go[[feature]], feature)) %>% 
 ug() %>% 
 mu(example = paste0(clean_name, "\n", cohortGo))

cohort_select <- "Pan-Cancer / Immunotherapy"

share <- 
next_step %>% 
 mu(focus = 
    (cohortGo == cohort_select & p_fdr_fisher_by < .1 & Odds == "Worse" & p_fdr_surv_by < .1 & !is.na(p_fdr_surv)), 
    Treatment = ifelse(focus, mechanism, "Other"))

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
  "hotspot_KRAS_G12D" = "KRAS G12D hotspot",
  "hotspot_KRAS_G12" = "KRAS G12(D/V/A/C/G) hotspot",
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

share %>% fi(grepl("effector", feature))

saveRDS(share, paste0(SHARE_DIR, "3_ready.rds"))
