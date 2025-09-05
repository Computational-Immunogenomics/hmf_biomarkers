source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

library(forcats)

META_DIR <- paste0(I_DIR, 'metadata/')

meta <- 
fread( paste0( META_DIR, "metadata_update_feb21_2025.csv")) %>% 
 se(sampleId, primaryTumorLocation, primaryTumorType, primaryTumorExtraDetails)

#cuppa_update <- 
#fread("/mnt/petasan_immunocomp/datasets/hartwig/cuppa/cuppa_update.csv") %>% 
# fi(rank == 1, clfName %in% c("COMBINED", "DNA_COMBINED")) 

cuppa <- 
fread("/mnt/petasan_immunocomp/datasets/hartwig/cuppa/cuppa_update_processed.csv") %>%
 mu(cuppaPrediction = ifelse(pred_class_combined_1 == "", pred_class_dna_1, pred_class_combined_1 ), 
    cuppaProb = ifelse(is.na(pred_prob_combined_1), pred_prob_dna_1, pred_prob_combined_1 )) %>% 
 se(sampleId, cuppaPrediction, cuppaProb) 

go <- 
meta %>% 
 lj(cuppa, by = "sampleId") %>% 
 mu(primaryTumorLocation = ifelse(primaryTumorLocation %in% c("NULL", "Unknown", "Other"), "Unspecified", primaryTumorLocation), 
    cuppaProb = ifelse(is.na(cuppaProb), 0, cuppaProb))

viewer <- function(i){
go %>% 
 fi(primaryTumorLocation %in% i) %>% 
 gb(primaryTumorLocation, primaryTumorType, primaryTumorExtraDetails, cuppaPrediction) %>% 
 su(ct = n()) %>% 
 ar(desc(ct))
}

breast_classifier <- function(primaryTumorExtraDetails){
 if(is.na(primaryTumorExtraDetails)){"Breast Unknown/Other"}
 else if( primaryTumorExtraDetails == "ER-negative/Her2-negative (triple negative)"){ "Triple Negative Breast Carcinoma"} 
 else if( primaryTumorExtraDetails == "ER-positive/Her2-negative" ){ "ER+/HER- Breast Carcinoma" } 
 else if( primaryTumorExtraDetails == "ER-positive/Her2-positive"){ "ER-/HER+ Breast Carcinoma"}      
 else if( primaryTumorExtraDetails == "ER-negative/Her2-positive"){ "ER+/HER+ Breast Carcinoma" } 
 else {"Unknown/Other Breast Carcinoma"}   
}
lung_classifier <- function(primaryTumorType, cuppaPrediction, cuppaProb){
 #print(primaryTumorType)
 if( primaryTumorType == "Non-small cell carcinoma (NSCLC), adenocarcinoma" ){ "Lung Adenocarcinoma" } 
 else if( primaryTumorType == "Non-small cell carcinoma (NSCLC), squamous cell carcinoma" ){ "Lung Squamous Cell Carcinoma" }   
 else if( primaryTumorType == "Non-small cell carcinoma (NSCLC), large cell neuroendocrine carcinoma" ){ "Lung Large Cell Carcinoma" }      
 else if( primaryTumorType == "Small cell carcinoma (SCLC)" ){ "Small Cell Lung Cancer" }   
 else if( grepl("Neuroendocrine tumor", primaryTumorType) ){ "Lung Neuroendocrine" }   
 else if( primaryTumorType == "Non-small cell carcinoma (NSCLC), not otherwise specified" & cuppaProb > .9){
  if(grepl("LUAD", cuppaPrediction)){ "Lung Adenocarcinoma" } 
  else if(grepl("LUSC", cuppaPrediction)){ "Lung Squamous Cell Carcinoma" } 
  else if(grepl("Small cell", cuppaPrediction)){ "Small Cell Lung Cancer" }
  else { "Lung Unknown/Other" } 
 }    
 else {"Lung Unknown/Other"}
}
skin_classifier <- function(primaryTumorType){
 primaryTumorType = tolower(primaryTumorType)
 if(grepl("melanoma", primaryTumorType)){ "Skin Melanoma"} 
 else if (grepl("basal cell carcinoma", primaryTumorType)){ "Skin Basal Cell Carcinoma"} 
 else if (grepl("squamous cell carcinoma", primaryTumorType)){ "Skin Squamous Cell Carcinoma"}    
 else if (grepl("merkel cell carcinoma", primaryTumorType)){ "Skin Merkel Cell Carcinoma"} 
 else {"Skin Other"}
}
soft_tissue_classifier <- function(primaryTumorType){
 primaryTumorType = tolower(primaryTumorType)  
 if(grepl("leiomyosarcoma", primaryTumorType)){ "Leiomyosarcoma"} 
 else if (grepl("liposarcoma", primaryTumorType)){ "Liposarcoma"} 
 else if (grepl("astrointestinal", primaryTumorType)){ "Gastrointestinal Stromal Tumor"}
 else if (grepl("undifferentiated", primaryTumorType)){"Soft tissue Undifferentiated"}
 else if (primaryTumorType == "myxofibrosarcoma"){"Myxofibrosarcoma"}
 else if (primaryTumorType == "angiosarcoma"){"Angiosarcoma"}
 else if (primaryTumorType == "synovial sarcoma"){"Synovial Sarcoma"}   
 else if (primaryTumorType == "solitary fibrous tumor"){"Soft tissue Solitary Fibrous Tumor"}      
 else {"Soft Tissue Sarcoma Unknown/Other"}   
}
pancreas_classifier <- function(primaryTumorType){
 if( primaryTumorType == "Neuroendocrine tumor (NET)" ){ "Pancreatic Neuroendocrine"} 
 else if(grepl("Adenocarcinoma", primaryTumorType) | grepl("Unknown", primaryTumorType)){ "Pancreatic Adenocarcinoma"}    
 else {"Pancreas Other"}
}
bladder_classifier <- function(primaryTumorType){
 primaryTumorType = tolower(primaryTumorType)
 if(grepl("urothelial carcinoma", primaryTumorType) | primaryTumorType == "unknown"){ "Urothelial Tract Bladder Adenocarcinoma"} 
 else {"Bladder Other"}
}
kidney_classifier <- function(primaryTumorType){
 if(grepl("Clear cell renal cell carcinoma", primaryTumorType)){ "Renal Clear Cell Carcinoma" } 
 #else if(grepl("Papillary renal cell carcinoma", primaryTumorType)){ "Renal non-clear cell carcinoma" }   
 else {"Renal (Unknown/Other)"}
}
colorectum_classifier <- function(primaryTumorType){
 if( primaryTumorType %in% c("Unknown", "Adenocarcinoma", "Mucinous adenocarcinoma") ){ "Colorectal Adenocarcinoma" }       
 else {"Colorectal Other"}
}
brain_classifier <- function(primaryTumorType){
 if( primaryTumorType == "Glioblastoma multiforme (GBM)" ){ "Glioblastoma Multiforme" } 
 else if( primaryTumorType == "Astrocytoma" ){ "Pilocityc Astrocytoma" } 
 else {"Brain/CNS Unknown/Other"}
}
cup_classifier <- function(cuppaPrediction, cuppaProb){
  if(is.na(cuppaPrediction)) { "CUPPA unknown"}
  else{ paste0("CUPPA predicted: ", cuppaPrediction) }
} 

cohort_map <- 
c("Prostate" = "Prostate adenocarcinoma", 
  "Ovary" = "Ovarian adenocarcinoma",
  "Biliary tract" = "Biliary tract adenocarcinoma", 
  "Liver" = "Hepatocellular carcinoma",
  "Pleura" = "Mesothelioma",
  "Ureter/renal pelvis" = "Urothelial tract non-bladder adenocarcinoma",
  "Urachus" = "Urothelial tract non-bladder adenocarcinoma",
  "Urethra" = "Urothelial tract non-bladder adenocarcinoma",
  "Esophagus/gastroesophageal junction" = "Gastroesophageal adenocarcinoma", 
  "Oral cavity/tongue" = "Head and Neck squamous cell carcinoma", 
  "Paranasal sinuses" = "Head and Neck squamous cell carcinoma",
  "Nasal cavity" = "Head and Neck squamous cell carcinoma",
  "Oropharynx" = "Head and Neck squamous cell carcinoma",
  "Nasopharynx" = "Head and Neck squamous cell carcinoma",
  "Hypopharynx" = "Head and Neck squamous cell carcinoma",
  "Larynx" = "Head and Neck squamous cell carcinoma", 
  "Gallbladder" = "Gallblader adenocarcinoma", 
  "Endometrium/uterus" = "Endometrium cancer", 
  "Thyroid gland" = "Thyroid adenocarcinoma", 
  "Vulva" = "Vulva/Vaginal squamous cell carcinoma",
  "Vagina" = "Vulva/Vaginal squamous cell carcinoma",
  "Penis" = "Penile cancer", 
  "Testis" = "Testicular cancer", 
  "Eye" = "Uveal Melanoma", 
  "Ampulla of Vater" = "Ampullary Cancer", 
  "Cervix" = "Cervix Carcinoma", 
  "Small bowel/intestine" = "Small bowel cancer", 
  "Anus" = "Anal Carcinoma", 
  "Salivary gland" = "Salivary gland carcinoma",
  "Bone" = "Bone Sarcoma", 
  "Appendix" = "Appendix Carcinoma", 
  "Thymus" = "Thymic Carcinoma",
  "Adrenal gland" = "Adrenal gland carcinoma",
  "Peritoneum" = "Mesothelioma",
  "Stomach" = "Stomach adenocarcinoma",
  "Lymphoid" = "Lymphoma",
  "Bone marrow/myeloid" = "Myeloid Cancer"
  )
cohort_mapper <- function(primaryTumorLocation){
  if(primaryTumorLocation %in% names(cohort_map)){str_to_title(cohort_map[[primaryTumorLocation]])}
  else { primaryTumorLocation }  
}

cohort_classifier <- function( primaryTumorLocation, primaryTumorType, primaryTumorExtraDetails, cuppaPrediction, cuppaProb ){
 if(is.na(primaryTumorLocation) | primaryTumorLocation == "NULL") {primaryTumorType}
 else if( primaryTumorLocation == "Breast" ){ breast_classifier(primaryTumorExtraDetails)}
 else if ( primaryTumorLocation == "Lung" ){ lung_classifier(primaryTumorType, cuppaPrediction, cuppaProb) }
 else if ( primaryTumorLocation == "Soft tissue" ) { soft_tissue_classifier(primaryTumorType) }
 else if ( primaryTumorLocation == "Skin" ) { skin_classifier(primaryTumorType) }
 else if ( primaryTumorLocation == "Bladder" ) { bladder_classifier(primaryTumorType) }   
 else if ( primaryTumorLocation == "Pancreas" ) { pancreas_classifier(primaryTumorType) } 
 else if ( primaryTumorLocation == "Kidney" ) { kidney_classifier(primaryTumorType) }    
 else if ( primaryTumorLocation %in% c("Colon", "Rectum") ) { colorectum_classifier(primaryTumorType) }    
 else if ( primaryTumorLocation == "Brain/central nervous system" ) { brain_classifier(primaryTumorType) }
 else if ( primaryTumorLocation %in% c("Cancer of unknown primary (CUP)") ) { cup_classifier(cuppaPrediction, cuppaProb) }
 else { cohort_mapper(primaryTumorLocation) }
}

acronyms <- 
c( 'Colorectal Adenocarcinoma' = "COREAD", 
	'ER+/HER- Breast Carcinoma' = "ER+/HER-BC", 
	'Lung Adenocarcinoma' = "LUAD", 
	'Prostate Adenocarcinoma' = "PRAD",
	'Skin Melanoma' = "SKCM",
	'Ovarian Adenocarcinoma' = "OV",
	'Gastroesophageal Adenocarcinoma' = "ESCA",
	'Pancreatic Adenocarcinoma' = "PAAD",
	'Triple Negative Breast Carcinoma' = "TNBC",
	'Urothelial Tract Bladder Adenocarcinoma' = "BLCA",
	'Soft Tissue Sarcoma Unknown/Other' = "SARCOT",
	'Glioblastoma Multiforme' = "GBM",
	'Unknown/Other Breast Carcinoma' = "BRCAUN",
	'Lung Unknown/Other' = "LUOT",
	'Biliary Tract Adenocarcinoma' = "CHOL",
	'ER-/HER+ Breast Carcinoma' = "ER-/HER+BC",
	'CUPPA unknown' = "CUP-NA",
	'Renal Clear Cell Carcinoma' = "RCCC",
	'Cervix Carcinoma' = "CESC",
	'Small Bowel Cancer' = "SBC",
	'Endometrium Cancer' = "EC",
	'Unspecified' = "NOS",
	'Gastrointestinal Stromal Tumor' = "GIST",
	'CUPPA predicted: Lung: Non-small cell: LUAD' = "CUP-LUAD",
	'Leiomyosarcoma' = "LEI",
	'Small Cell Lung Cancer' = "SCLC",
	'Renal (Unknown/Other)' = "RCOT",
	'Urothelial Tract Non-Bladder Adenocarcinoma' = "UTNBC",
	'Stomach Adenocarcinoma' = "STAD",
	'Lung Neuroendocrine' = "LNET",
	'Mesothelioma' = "MESO",
	'Head And Neck Squamous Cell Carcinoma' = "HNSC",
	'Hepatocellular Carcinoma' = "HCC",
	'CUPPA predicted: HPB: Bile duct/Gallbladder' = "CUP-HPB",
	'Liposarcoma' = "LIPO",
	'Lung Squamous Cell Carcinoma' = "LUSC",
	'Soft tissue Undifferentiated' = "SARCUN",
	'Anal Carcinoma' = "ANCA",
	'CUPPA predicted: Colorectum/Small intestine/Appendix' = "CUP-CSA",
	'ER+/HER+ Breast Carcinoma' = "ER+HER+BC",
	'CUPPA predicted: HPB: Pancreas' = "CUP-HPB",
	'Brain/CNS Unknown/Other' = "BROT",
	'Pancreatic Neuroendocrine' = "PANET",
	'Colorectal Other' = "COOT",
	'Thyroid Adenocarcinoma' = "THCA",
	'CUPPA predicted: Esophagus/Stomach' = "CUP-ES",
	'Lymphoma' = "LYMPH",
	'Pilocityc Astrocytoma' = "PIA",
	'Salivary Gland Carcinoma' = "SGCA",
	'Gallblader Adenocarcinoma' = "GACA",
	'Skin Basal Cell Carcinoma' = "SKBC",
	'Vulva/Vaginal Squamous Cell Carcinoma' = "VVSC",
	'Appendix Carcinoma' = "APCA",
	'CUPPA predicted: Anogenital' = "CUP-ANG",
	'Bone Sarcoma' = "OS",
	'CUPPA predicted: Bone/Soft tissue: Other' = "CUP-BST",
	'CUPPA predicted: Gynecologic: Ovary/Fallopian tube' = "CUP-GYN",
	'Thymic Carcinoma' = "TYCA",
	'Lung Large Cell Carcinoma' = "LCLC",
	'CUPPA predicted: Skin: Melanoma' = "CUP-SKCM",
	'Penile Cancer' = "PNCA",
	'CUPPA predicted: NET: Pancreas' = "CUP-PANET",
	'Skin Other' = "SKOT",
	'Angiosarcoma' = "ANGS",
	'Bladder Other' = "BLOT",
	'Skin Squamous Cell Carcinoma' = "SKSC",
	'Soft tissue Solitary Fibrous Tumor' = "SARCFT",
	'Lung Squamous Cell Carcinoma' = "LUSC",
	'Synovial Sarcoma' = "SYN",
	'CUPPA predicted: Urothelial tract' = "CUP-UT",
	'Adrenal Gland Carcinoma' = "AACA",
	'CUPPA predicted: Head and neck: Other' = "CUP-HNSC",
	'CUPPA predicted: Kidney: Other' = "CUP-RCOT",
	'CUPPA predicted: Lung: Small cell' = "CUP-SCLC",
	'Myxofibrosarcoma' = "MYX",
	'Pancreas Other' = "PAOT",
	'CUPPA predicted: Breast: Other' = "CUP-BRCAUN",
	'CUPPA predicted: Skin: Other' = "CUP-SKOT",
	'Skin Merkel Cell Carcinoma' = "SKMC",
	'CUPPA predicted: Prostate' = "CUP-PRAD",
	'Testicular Cancer' = "TTCA",
	'CUPPA predicted: HPB: Liver' = "CUP-HPB",
	'CUPPA predicted: NET: Colorectum/Small intestine' = "CUP-CSNET",
	'CUPPA predicted: NET: Lung' = "CUP-LNET",
	'CUPPA predicted: Breast: Triple negative' = "CUP-TNBC",
	'CUPPA predicted: Gynecologic: Endometrium' = "CUP-GYNEND",
	'CUPPA predicted: Lung: Non-small cell: LUSC' = "CUP-LUSC",
	'Uveal Melanoma' = "UVM",
	'CUPPA predicted: Mesothelium' = "CUP-MESO",
	'Ampullary Cancer' = "AC",
	'CUPPA predicted: Bone/Soft tissue: Leiomyosarcoma' = "CUP-LEI",
	'CUPPA predicted: CNS: Glioma' = "CUP-GLI",
	'CUPPA predicted: Head and neck: Adenoid cystic' = "CUP-HNAC",
	'CUPPA predicted: Bone/Soft tissue: Undiff. sarcoma' = "CUP-SARCUN",
	'CUPPA predicted: Head and neck: Salivary gland' = "CUP-HNSG",
	'CUPPA predicted: Lymphoid tissue' = "CUP-LYMPH",
	'Myeloid Cancer' = "MC",
	'CUPPA predicted: CNS: Medulloblastoma' = "CUP-CNSMED",
	'CUPPA predicted: Kidney: Chromophobe' = "CUP-RCCC",
	'CUPPA predicted: Myeloid: Myeloproliferative neoplasm' = "CUP-MYL",
	'CUPPA predicted: Thyroid gland' = "CUP-THCA")

cohorts_ready <- 
go %>% 
 rw() %>% 
 mu(cohort = cohort_classifier(primaryTumorLocation, primaryTumorType, primaryTumorExtraDetails, cuppaPrediction, cuppaProb), 
    acronym = acronyms[[cohort]]) %>% 
 ug() %>% 
 unique() %>% 
 drop_na(cohort)

summaries <- 
cohorts_ready %>% 
 gb(primaryTumorLocation, cohort, acronym) %>% 
 su(ct = n()) %>% 
 gb(primaryTumorLocation) %>% mu(tot = sum(ct)) %>% ar(desc(tot), desc(ct))

fwrite(summaries, "summary.csv")

fwrite(cohorts_ready, paste0(META_DIR, "cohorts/cohorts_ready.csv"))
