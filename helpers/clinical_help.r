options(repr.matrix.max.cols=100, repr.matrix.max.rows=100)

date_diff <- function(d1, d2) {
  if( is.na(d1) || is.na(d2) || d1 == "NULL" || d2 == "NULL"){ NA } 
  else {as.numeric( difftime(d2, d1, units = "days"))}
}

recist_name_map <- 
list(
  'CR' = 'CR',
  'ICR' = 'CR',  
  'PR' = 'PR',
  'IPR' = 'PR', 
  'CLINICAL_BENEFIT' = 'PR',  
  'SD' = 'SD',
  'ISD' = 'SD', 
  'NON CR/NON PD' = 'SD',  
  'NON-CR/NON-PD' = 'SD',
  'NON ICR/NON IPD' = 'SD',  
  'PD' = 'PD',
  'IPD' = 'PD',
  'IUPD' = 'PD',  
  'CLINICAL PROGRESSION' = 'PD', 
  'STOP_TREATMENT;CLINICAL_DETERIORATION' = 'PD',
  'STOP_TREATMENT;RADIOLOGIC_PROGRESSION_NON_RECIST_LESIONS' = 'PD',
  'STOP_TREATMENT;VISIBLE_EVIDENT_PROGRESSION' = 'PD'
)

derive_response <- function(i){
  i <- toupper(i)
  if( i %in% names(recist_name_map) ){ recist_name_map[[i]]} 
  else { NA }
}

go_dcb <- function( i ){
  if( grepl("CR",i) || grepl("PR",i) || grepl("SD_durable", i) ){ 1 } 
  else if ( grepl("SD",i) || grepl("PD",i) ) { 0 }
  else { NA }
}

go_bor <- function( i ){
  if( grepl("CR",i) || grepl("PR",i)){ 1 } 
  else if ( grepl("SD",i) || grepl("PD",i) ) { 0 }
  else { NA }
}

min2 <- function(i, j){
 if(is.na(i) & is.na(j)){ NA } 
 else { min(i,j, na.rm = TRUE)}
}

nicer_treatment_format <- function(i) {
  paste0(sort(unique(strsplit(i, "/")[[1]])), collapse = "/")
}

yes_no_indicator <- function(i) {
 j <- tolower(i)
 if( j == "yes"){1}
 else if( j == "no"){0}
 else {NA}
}

get_year <- function(i) {
  as.numeric(strsplit(i, "-")[[1]][1])
}

sex <- function(gender){
  if(!is.na(gender)){
    if( tolower(gender) == "female"){ 1 }
    else { 0 }
  } else { NA }
}

age <- function( birthYear, biopsyDate ){
  if(!is.na(birthYear) & !is.na(biopsyDate)){
    biopsyYear <- as.numeric(strsplit(biopsyDate, "-")[[1]][1])
    age <- biopsyYear - as.numeric(birthYear)
  } else { NA }
}

responder <- function(v){
  if( "CR" %in% v || "PR" %in% v ){ 1 }
  else{ 0 }
}

trt_indicator <- function (v, t) {
  apply(data.frame(lapply(v, function(i) grepl(i, t))), 1, sum)
}
