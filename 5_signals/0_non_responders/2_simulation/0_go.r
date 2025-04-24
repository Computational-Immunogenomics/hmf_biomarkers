source(paste0(dirname(dirname(dirname(getwd()))),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

go <- function(n = 100, prevalence = .5, p_base = .4, p_event = .05) {

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
  non_responders <- n - sum(ready$response)  
    
  if( (events > 0) && (responders > 0)){  
    fisher <- fisher.test(table(ready))
    df( n = n, 
        events = events, 
        non_events = non_events, 
        responders = responders, 
        prevalence = prevalence, 
        p_base = p_base, 
        p_event = p_event,
        p_fisher = fisher$p.value, 
        or = fisher$estimate, 
        ci.low = fisher$conf.int[1], 
        ci.high = fisher$conf.int[2])
  } else {
    df( n = n, events = events, non_events = non_events, responders = NA, 
        prevalence = prevalence, p_base = p_base, p_event = p_event, p_fisher = 1, or = NA, ci.low = NA, ci.high = NA)
  }
}

set.seed(62220)
nsim <- 1000
ns <- c(20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 1000)
prevalence <- c(.01, .05, .1, .3, .5, .8)
p_base <- c(.4, .3, .2, .1)
#p_base <- c(.4)
p_event <- c(.4, .2, .1, .01, 0)

oo <- data.frame()
system.time(    
for( z in seq(nsim)){
 print(z); flush.console();
 for( i in ns){
  for( j in p_base ){  
   for ( k in p_event ){
    for( l in prevalence ) {   
     tmp <- tryCatch({go(n = i, prevalence = l, p_base = j, p_event = k)}, error = function(e) {return(NA)})
     if(is.data.frame(tmp)) oo <- rbind(oo, tmp)
}}}}})

fwrite(oo, "sim_go_new.csv")
