source(paste0(dirname(getwd()),'/map.r'))

source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "clinical_help.r"))

META_DIR <- paste0(I_DIR, 'metadata/')

arranger <- function(i) paste(unlist(list(sort(strsplit(unique(i), " ## ")[[1]]))), collapse = " ## ")

META_DIR

meta <- fread( paste0( META_DIR, "metadata_update_feb21_2025.csv"))
response <- fread( paste0( META_DIR, "treatment_responses.tsv"))
pre_biopsy_drugs <- fread( paste0( META_DIR, "pre_biopsy_drugs.tsv"))
post_biopsy_drugs <- fread( paste0( META_DIR, "post_biopsy_drugs.tsv"))
new_biopsies <- fread(paste0(META_DIR, "sophie_share/dr347_biopsyInfo.csv"))

head(meta)

pretreatments <- 
pre_biopsy_drugs %>% 
 gb(patientIdentifier) %>% 
 su( preTreatmentName = arranger(paste0(name, collapse = " ## ")), 
     preTreatmentType = arranger(paste0(type, collapse = " ## ")),
     preTreatmentMechanism = arranger(paste0(mechanism, collapse = " ## ")),
     preTreatmentLines = n_distinct(startDate))

post_biopsy_treatments <- 
post_biopsy_drugs %>% 
 mu(treatment = paste0(sampleId,"##", startDate),
    startDate = as.character(startDate)) %>% 
 ug() %>% 
 gb(treatment) %>% 
 su(patientIdentifier = unique(patientIdentifier),
    sampleId = unique(sampleId), 
    treatmentStartDate = min(startDate), 
    treatmentEndDate = max(endDate), 
    treatmentName = arranger(paste0(name, collapse = " ## ")), 
    treatmentType = arranger(paste0(type, collapse = " ## ")), 
    treatmentMechanism = arranger(paste0(mechanism, collapse = " ## "))) %>% 
 gb(patientIdentifier) %>% 
 mu(postBiopsyTreatmentLine = rank(treatmentStartDate)) %>% 
 ug() 

outcomes <- 
response %>% 
 mu(treatment = paste0(sampleId,"##", as.character(startDate)), 
    responseDate = as.character(responseDate)) %>% 
 gb(patientIdentifier) %>% 
 mu(postBiopsyTreatmentLine = dense_rank(startDate)) %>% 
 ug()

treatments_and_outcomes <- 
pretreatments %>% 
 full_join(post_biopsy_treatments, by = "patientIdentifier") %>% 
 full_join(outcomes %>% se(treatment, responseDate, response), by = "treatment") %>% 
 relocate(treatment, patientIdentifier, sampleId)

ageify <- function( birthYear, biopsyDate ) {
  if(birthYear == "NULL" || is.na(biopsyDate)){ NA } 
  else { as.numeric(substr(biopsyDate, 0, 4)) - as.numeric(birthYear) }
}
genderify <- function( gender ) {
  if( gender == "NULL" ){ NA } 
  else if (gender == "female") { 1 }
  else { 0 }  
}

meta_ready <- 
meta %>% 
 rename(biopsyTypeOld = biopsyType, biopsyDateOld = biopsyDate) %>% 
 lj( new_biopsies %>% 
     tm(sampleId, biopsyStructure, biopsyType, biopsyDate = format(as.Date(biopsyDate, format = "%d-%m-%Y"), "%Y-%m-%d")), 
  by = "sampleId") %>% 
 rw() %>% 
 mu( age = ageify(birthYear, biopsyDate), sex = genderify( clinicalGender ) )

metadata_dates <- meta_ready %>% se( patientId, sampleId, sampleArrivalDate, biopsyDate, deathDate) 

date_diff <- function (d1, d2) {
 if (is.na(d1) || is.na(d2) || tolower(d1) == "null" || tolower(d2) == "null") { NA } 
 else { as.numeric(difftime(d2, d1, units = "days"))}
}

together <- 
metadata_dates %>% 
 lj(treatments_and_outcomes, by = "sampleId") %>% 
 rw() %>% 
 mu(
  raw_response = response,  
  response = derive_response(response),                             ### applying the recist_name_map
  os_event = ifelse(deathDate  != "NULL", 1, 0),                    ### do we have a death date
  pfs_event = ifelse(os_event == 1 || response == "PD", 1, 0),      ### death or progression
  days_to_treatment = date_diff(biopsyDate, treatmentStartDate), 
  days_to_treatment_end = date_diff( as.character(treatmentStartDate), as.character(treatmentEndDate)),
  days_to_response = date_diff( treatmentStartDate, responseDate ),
  days_to_death = date_diff( treatmentStartDate, deathDate ) , 
  days_to_last_measured = max( days_to_response, days_to_treatment_end, days_to_death, na.rm = TRUE ),
  days_to_progression = ifelse( response == "PD", days_to_response, NA ),
  days_to_pfs = ifelse(pfs_event == 1, min2(days_to_progression, days_to_death), NA), 
  response_dcb = ifelse(days_to_response >= 183 & response == "SD", "SD_durable", response)
  ) %>% 
 ug() %>% 
 mu(days_to_last_measured = ifelse(days_to_last_measured == "-Inf", NA, days_to_last_measured))

clinical_outcomes <- 
together %>% 
 gb(patientId, sampleId, treatmentId = treatment) %>%              ### summarise at treatment Id levels
 su(
  derived_preTreatmentName = unique(preTreatmentName), 
  derived_preTreatmentType = unique(preTreatmentType), 
  derived_preTreatmentMechanism = unique(preTreatmentMechanism), 
  derived_preTreatmentLines	 = unique(preTreatmentLines), 
  derived_treatmentName = unique(treatmentName),
  derived_treatmentType = unique(treatmentType), 
  derived_treatmentMechanism = unique(treatmentMechanism),
  rawResponses = paste0(unique(raw_response), collapse = ","), 
  responses = paste0(unique(response_dcb), collapse = ","),  
  completeResponse = ifelse(grepl("CR", responses),1,0),
  completeResponse = ifelse(is.na(responses), NA, completeResponse),    
  bestOverallResponse = go_bor(responses), 
  durableClinicalBenefit = go_dcb(responses), 
  pfsEvent = max(pfs_event, na.rm = TRUE),
  daysToPfsEvent = ifelse(pfsEvent==1, min(days_to_pfs, na.rm=TRUE), max(days_to_last_measured, na.rm=TRUE)),
  osEvent = max(os_event, na.rm = TRUE), 
  daysToOsEvent = max(days_to_last_measured, na.rm = TRUE),
  postInitialBiopsyTreatmentLine = mean( postBiopsyTreatmentLine, na.rm = TRUE), 
  daysBiopsyToTreatment = mean(days_to_treatment, na.rm = TRUE)) %>% 
 ug() 

clin_base <- 
meta_ready %>% 
 se(-contains("date")) %>% 
 full_join(clinical_outcomes, by = c("sampleId", "patientId")) %>% 
 mu(across(everything(), ~ replace(., . %in% c(-Inf, NaN, NULL, "NULL"), NA)))

fwrite( clin_base, paste0(TMP_DIR, "clinical.csv") )

clin_ready <- 
clin_base %>% 
 se(patientId, sampleId, 
    primaryTumorLocation, primaryTumorType, 
    age, sex, 
    contains("derived"), contains("therapy"), hasSystemicPreTreatment, biopsyStructure, biopsyType, 
    bestOverallResponse, durableClinicalBenefit, 
    pfsEvent, daysToPfsEvent, 
    osEvent, daysToOsEvent,
    postInitialBiopsyTreatmentLine, daysBiopsyToTreatment) %>%
    mu(hasSystemicPreTreatment = as.numeric(hasSystemicPreTreatment == "Yes"), 
       hasRadiotherapyPreTreatment	 = as.numeric(hasRadiotherapyPreTreatment == "Yes"),
       radiotherapyGiven	 = as.numeric(radiotherapyGiven == "Yes")) %>% 
 gb(sampleId) %>% mu(rk = row_number(daysBiopsyToTreatment)) %>% ug() %>% 
 mu(rk = ifelse(is.na(rk), 1, rk)) %>% fi(rk == 1) %>% se(-rk) %>% 
 rename_with( ~paste0("clin_", .x), 
    .cols = c("age", "sex", "hasRadiotherapyPreTreatment", "radiotherapyGiven", "hasSystemicPreTreatment", 
              "postInitialBiopsyTreatmentLine", "daysBiopsyToTreatment"))

mechanisms <- pre_biopsy_drugs %>% gb(mechanism) %>% su(ct = n()) %>% ar(desc(ct)) %>% fi(ct > 30) %>% pu(mechanism)

for(i in mechanisms) {
 clin_ready[,paste0("clin_pretreatment_contains_", gsub("[^a-zA-Z]", "", i))] <- 
 clin_ready %>% mu( tmp = as.numeric(grepl(i, derived_preTreatmentMechanism))) %>% pu(tmp)
}

biopsy_sites <- 
clin_ready %>% 
 gb(biopsyStructure) %>% 
 su(ct = n()) %>% drop_na() %>% fi(biopsyStructure != "Not specified", ct > 30) %>% 
 pu(biopsyStructure)

for(i in biopsy_sites) {
 clin_ready[,paste0("clin_biopsy_site_", gsub("[^a-zA-Z]", "", i))] <- 
 clin_ready %>% mu( tmp = as.numeric(biopsyStructure == i)) %>% pu(tmp)
}

clin_ready <- 
clin_ready %>% 
 mu( clin_biopsy_type = case_when(
      biopsyType == "Not specified" ~ NA,
      biopsyType == "Metastasis" ~ 1,
      biopsyType == "Primary" ~ 0))

fwrite( clin_ready, paste0(READY_DIR, "clinical_ready.csv") )
