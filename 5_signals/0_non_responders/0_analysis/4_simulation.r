source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

go <- function(n = 100, prevalence = .5, p_base = .4, p_event = .05, z = 1) {

  ### Simulate events 
  X <- rbinom(n, 1, prevalence)
  events <- sum(X)
  non_events <- n - events

  ### Simulate response   
  if( events > 0 ){
    Y_event <- df( "event" = 1, response = rbinom( events, 1, p_event))
    Y_non_event <- df( "event" = 0, response = rbinom( non_events, 1, p_base))
    ready <- rbind(Y_event, Y_non_event) 
  } else {
    ready<- df( "event" = 0, response = rbinom( non_events, 1, p_base))
  }
  responders <- sum(ready$response)
  responders_event <- ready %>% mu(event_and_response = ifelse(event + response == 2, 1, 0)) %>% su(a = sum(event_and_response)) %>% pu(a)
  non_responders <- n - sum(ready$response)  
    
  if( (events > 0) && (responders > 0)){  
    fisher <- fisher.test(table(ready))
    df( n = n, 
        events = events, 
        responders_event = responders_event, 
        non_events = non_events, 
        responders = responders, 
        prevalence = prevalence, 
        p_base = p_base, 
        p_event = p_event,
        p_fisher = fisher$p.value, 
        or = fisher$estimate, 
        ci.low = fisher$conf.int[1], 
        ci.high = fisher$conf.int[2], 
        z = z)
  } else {
    df( n = n, events = events, responders_event = NA, non_events = non_events, responders = NA, 
        prevalence = prevalence, p_base = p_base, p_event = p_event, p_fisher = 1, or = NA, ci.low = NA, ci.high = NA, z = z)
  }
}

set.seed(62220)
nsim <- 1000
ns <- c(30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 1000)
prevalence <- c(.01, .05, .1, .5)
p_base <- c(.4, .2)
p_event <- c(.4, .2, .1, .01, 0)

oo <- data.frame()
system.time(    
for( z in seq(nsim)){
 print(z); flush.console();
 for( i in ns){
  for( j in p_base ){  
   for ( k in p_event ){
    for( l in prevalence ) {   
     if( k <= j){
     tmp <- tryCatch({go(n = i, prevalence = l, p_base = j, p_event = k, z = z)}, error = function(e) {return(NA)})
     rownames(tmp) <- NULL
     if(is.data.frame(tmp)) oo <- rbind(oo, tmp)
}}}}}})

go_binom_test <- function( n, x, p = .02) {
   if(is.na(x)){1}
   else if (n == 0) {1}
   else{ binom.test(x, n, p, alternative = "less")$p.value } 
}

go <- 
oo %>% 
 rw() %>% 
 mu(expected_events = n*prevalence, 
    pval_under01 = go_binom_test(events, responders_event, .01), 
    pval_under2 = go_binom_test(events, responders_event, .02), 
    pval_under5 = go_binom_test(events, responders_event, .05),
    pval_under10 = go_binom_test(events, responders_event, .1)) %>% 
 ug()

fwrite(oo, paste0(SHARE_DIR, "5_simulation_results.csv"))
