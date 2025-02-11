source(paste0(dirname(getwd()),'/map.r'))

source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "clinical_help.r"))

META_DIR <- paste0(I_DIR, 'metadata/')

meta <- fread( paste0( META_DIR, "metadata_update.csv"))
response <- fread( paste0( META_DIR, "treatment_responses.tsv"))
pre_biopsy_drugs <- fread( paste0( META_DIR, "pre_biopsy_drugs.tsv"))
post_biopsy_drugs <- fread( paste0( META_DIR, "post_biopsy_drugs.tsv"))
new_biopsies <- fread(paste0(META_DIR, "sophie_share/dr347_biopsyInfo.csv"))

pretreatments <- 
pre_biopsy_drugs %>% 
 gb(patientIdentifier) %>% 
 su( preTreatmentName = paste0(unique(name), collapse = "/"), 
     preTreatmentType = paste0(unique(type), collapse = "/"),
     preTreatmentMechanism = paste0(unique(mechanism), collapse = "/"),
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
    treatmentName = paste0(unique(name), collapse = "/"), 
    treatmentType = paste0(unique(type), collapse = "/"), 
    treatmentMechanism = paste0(unique(mechanism), collapse = "/")) %>% 
 gb(patientIdentifier) %>% 
 mu(postBiopsyTreatmentLine = rank(treatmentStartDate)) %>% 
 ug()

outcomes <- 
response %>% 
 mu(treatment = paste0(sampleId,"##", as.character(startDate)), 
    responseDate = as.character(responseDate)) %>% 
 gb(patientIdentifier) %>% 
 mu(postBiopsyTreatMentLine = dense_rank(startDate)) %>% 
 ug()

treatments_and_outcomes <- 
pretreatments %>% 
 full_join(post_biopsy_treatments, by = "patientIdentifier") %>% 
 full_join(outcomes %>% se(treatment, responseDate, response), by = "treatment") %>% 
 relocate(treatment, patientIdentifier, sampleId)

meta_ready <- 
meta %>% 
 rename(biopsyTypeOld = biopsyType, biopsyDateOld = biopsyDate) %>% 
 lj( new_biopsies %>% 
     tm(sampleId, biopsyStructure, biopsyType, biopsyDate = format(as.Date(biopsyDate, format = "%d-%m-%Y"), "%Y-%m-%d")), 
  by = "sampleId")

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

clin_ready <- 
meta_ready %>% 
 full_join(clinical_outcomes, by = c("sampleId", "patientId")) %>% 
 mu(across(everything(), ~ replace(., . %in% c(-Inf, NaN, NULL, "NULL"), NA)))

anti_pd_trts <- c("Atezolizumab", "Avelumab", "Durvalumab", "Nivolumab", "Pembrolizumab")

clin_ready$clin_anti_pd_treated <- ifelse(trt_indicator(anti_pd_trts, clin_ready$clin_treatment) > 0, 1, 0)

fwrite( clin_ready, paste0(READY_DIR, "clinical_ready.csv") )

sarcoma_samples <- 
read.csv("/mnt/bioinfnas2/immunocomp/shared_reference_data/cohorts/sarcoma_samples.csv") %>% 
 rw() %>% 
 mu(samples = strsplit(sampleId, ",")[[1]][1]) %>% 
 pu(samples)

fwrite( clin_ready %>% fi(sampleId %in% sarcoma_samples), paste0(READY_DIR, "clinical_sarcoma_ready.csv") )

#paste0(READY_DIR, "clinical_sarcoma_ready.csv")

#clin_ready %>% fi(sampleId %in% sarcoma_samples) %>% fi(!is.na(responseMeasured))
