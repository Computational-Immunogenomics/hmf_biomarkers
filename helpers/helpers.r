library(survival)
get_model_formula <- function( y, x, full = TRUE){
  if( full ){
    as.formula(paste(y, x, sep = "~"))
  } else {
    tmp <- strsplit(x,"+",fixed = T)[[1]]
    x_reduced <- ifelse(length(tmp) == 1, "1", paste(tmp[-1], collapse='+' ))
    as.formula(paste(y, x_reduced, sep = "~"))
  }  
}
get_model <- function(y, x, covariate = "", data = "df", model = "bor", full_model = TRUE ){
  formula <- get_model_formula(y, paste0(x, covariate), full = full_model)
  if (model == "lm"){ do.call("glm", list(formula = formula, data = as.name(data), family = "gaussian" ))}
  else if (model == "bor"){ do.call("glm", list(formula = formula, data = as.name(data), family = "binomial" ))}
  else if ( model == "coxph"){ do.call("coxph",list(formula = formula,data = as.name(data)))}
}
coefs_extract <- function(model, type = "glm"){
  o <- summary(model)
  if (type == "coxph") {coefs <- as.numeric(o$coefficients[1,c(1,3,5)])}
  else { coefs <- as.numeric(o$coefficients[2,c(1,2,4)]) }
  list("est" = coefs[1], "se" = coefs[2], "pval" = coefs[3])
}
get_lrt_pval <- function( full, red, type = "coxph", covariate = "" ){
  if( type == "coxph" && covariate == ""){
      list( "lrt_pval" = NA)
  } else {
    if( type == "coxph" ){  df = as.numeric(summary(full)$logtest[2]) - as.numeric(summary(red)$logtest[2])  }
    else { df <- summary(red)$df.residual  - summary(full)$df.residual }
    test_stat <- 2*(as.numeric(logLik(full)) - as.numeric(logLik(red)))
    list("lrt_pval" = 1-pchisq(test_stat, df ) )
  }
}
get_stats <- function( y, x, covariate = "", data = "df", model = "bor"){
  full <- get_model(y = y, x = x, covariate = covariate, data = data, model = model)
  red <- get_model(y = y, x = x, covariate = covariate, data = data, model = model, full_model = FALSE)
  x_type <- strsplit(x, "_")[[1]][1]
  settings <- list("y" = y, "x" = x, "covariate" = covariate, "type" = x_type, "data" = data, "model" = model)
  coefs <- coefs_extract(full, type = model)
  lrt_pval <- get_lrt_pval(full, red, type = model, covariate = covariate)
  data.frame(append(append(settings, coefs), lrt_pval))
}
get_stats2 <- function( y, x, covariate = "", data = "df", model = "bor" ){
  out <- tryCatch( get_stats( y, x, covariate, data, model), error = function(e) NULL)
  return(out)
}
scanner <- function (y = "Surv(Y_os_days, Y_os_event)", features, covariates, df = "df", mod = "coxph") {
    out <- data.frame()
    for (f in features) {
        if (grepl("rna_", f)) {
            tmp <- get_stats2(y = y, x = f, covariate = paste0("+purity", covariates), data = df, model = mod)
        }
        else {
            tmp <- get_stats2(y = y, x = f, covariate = covariates, data = df, model = mod)
        }
        if (is.data.frame(tmp)) out <- rbind(out, tmp)
    }
    out
}
mod_map <- list(
 "bestOverallResponse" = "bor",
 "durableClinicalBenefit" = "bor",
 "Surv(daysToOsEvent, osEvent)" = "coxph",
 "Surv(daysToPfsEvent, pfsEvent)" = "coxph",
 "nr_bor" = "bor", 
 "nr_dcb" = "bor"
)