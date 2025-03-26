source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

go <- function(n = 100, beta = -1, pred = "bin", pp = .5, b0 = -1, cov = TRUE) {
  # Simulate predictor and response
  if(pred == "bin") { X <- rbinom(n, 1, pp)}
  else { X <- rnorm(n) }
  
  if(cov) { 
    X1 <- rbinom(n, 1, .3); X2 <- rbinom(n, 1, .7); X3 <- rnorm(n)
    cov_effect <- X1*1 + X2*(-1) + X3*.5
    log_odds <- b0 + beta*X + cov_effect
  } else { 
    log_odds <- b0 + beta*X 
  }
  
  p <- 1/(1 + exp(-log_odds)) ### Probability
  Y <- rbinom(n, 1, p) # Binary response
  
  # Fit logistic regression, with or without covariates
  if(cov){ model <- glm(Y ~ X + X1 + X2 + X3, family = binomial)} 
  else { model <- glm(Y ~ X, family = binomial)}
    
  # Extract the output for coefficient X 
  oo <- data.frame(t(summary(model)$coefficients[2,])) # p-value for X
    
  # Metadata   
  responders = sum(Y)

  if(pred == "bin"){
    p_event = mean(1/(1 + exp(-(b0+beta))))  #mean(p[which(X == 1)], na.rm = TRUE)
    p_non_event = mean(1/(1 + exp(-b0)))
  } else {
    p_event <- 1 / (1 + exp(-b0 + beta*1.96))
    p_non_event <- 1 / (1 + exp(-b0 - beta*1.96))  
  }
    
  p_overall = mean(p, na.rm = TRUE) 
  events = sum(X == 1)
  responders_event = sum(Y[which(X == 1)], na.rm = TRUE)
  p_simple <- pbinom(responders_event, events, p_overall)
  p_fisher <- fisher.test(table(data.frame(X, Y)))$p.value
    
  oo %>% 
   mu( n = n, 
       beta = beta, 
       pred = pred, 
       pp = pp, 
       b0 = b0, 
       responders = responders, 
       p_event = p_event,
       p_non_event = p_non_event,
       p_overall = p_overall, 
       events = events,
       responders_event = responders_event,
       p_simple = p_simple,
       p_fisher = p_fisher,
       cov = cov)
}

go_sim <- function(n = 100, beta = -1, pred = "bin", pp = .5, b0 = -1, cov = TRUE) {
  tryCatch({ 
    oo <- go(n, beta, pred, pp, b0, cov)
    return(oo)
  }, error = function(e) {
    return(NULL) # Return NULL or handle as needed
  })
}

set.seed(62220)

nsim <- 1000
ns <- c(25, 50, 100, 200, 500)
betas <- -c(0, 1.3, 2.6, 3.9)
preds <- c("bin")
pps <- c(.3, .1, .5, .05)
b0s = c(-.5)
covs = c(FALSE, TRUE)

oo <- data.frame()
system.time(    
for( z in seq(nsim)){
 print(z); flush.console();
 for( i in ns){
  for( j in betas ){  
   for ( k in preds){
    for( l in pps){
     for( q in b0s){
      for (c in covs){
       if( k == "norm"){oo <- rbind(oo, go_sim(n = i, beta = j/2, pred = k, pp = l, b0 = q, cov = c)) } 
       else { oo <- rbind(oo, go_sim(n = i, beta = j, pred = k, pp = l, b0 = q, cov = c))}
}}}}}}})

fwrite(oo, "sim_go.csv")
