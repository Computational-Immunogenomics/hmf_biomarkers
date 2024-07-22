sex <- function(gender){
  if(!is.na(gender)){
    if( tolower(gender) == "female"){ 1 }
    else { 0 }
  } else{ NA }
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
