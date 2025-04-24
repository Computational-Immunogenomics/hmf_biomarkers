require(dplyr)

min_patients <- 40
min_response <- 15
min_events <- 3
combination_threshold <- .01

ra_fisher <- function(a,b,c,d){
 fisher.test(matrix(c(a,b,c,d), ncol = 2))
}

ra_formatter_and_test <- function(df = "DF"){
 ### require df to be dataframe with cohort, non_response indicator, and columns of binary features
 df %>%
 ga(feature, event, -cohortGo, -non_response) %>%
 drop_na(event) %>%
 mu(non_response = ifelse(non_response == 1, "nr", "r"), event = ifelse(event == 1, "e", "ne")) %>%
 gb(cohortGo, feature, non_response, event) %>%
 su(tot = n(), .groups = "drop") %>%
 complete(cohortGo, feature, non_response, event, fill = list(tot = 0)) %>%
 pivot_wider(names_from = c(event, non_response),  values_from = tot) %>%
 mu(across(everything(), ~replace_na(., 0)),
    events = e_r + e_nr, no_events = ne_r + ne_nr, responders = e_r + ne_r, non_responders = e_nr + ne_nr, total_patients = events + no_events) %>%
 fi( events > min_events, no_events > min_events ) %>%
 mu(direction = ifelse( e_nr/events > non_responders/total_patients, "Non-Response", "Response")) %>%
 rw() %>%
  mu( results = list({
    oo = ra_fisher(`ne_nr`, `ne_r`, `e_nr`, `e_r`)
    tibble(fisher_pval = oo$p.value, or = oo$estimate, ci_low = oo$conf.int[1], ci_high = oo$conf.int[2])
  })) %>%
  unnest_wider(results)
}
