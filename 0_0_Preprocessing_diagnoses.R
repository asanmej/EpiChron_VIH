# Author: Alejandro Santos Mejías
# Date last update: 2025-04-25
# Input: Original bdu and diagnosis data
# Output: First draft bdu and transformed diagnosis from ICPC to CCS
# Motivational comment: Even from scratches something can be made!!!

source("99_paths_and_packages.R", encoding = "UTF-8")

##################################
##### Preprocessing Diagnoses ####
##################################

################################################################################
# Load dictionaries:
ccs <- fread("ccs_r/ccs_icd9_dic_vih.csv", colClasses = "character", encoding = "UTF-8")
cci <- fread("ccs_r/cci2015.csv", quote = "", encoding = "UTF-8")
cohort <- fread("control/controles.csv", encoding = "UTF-8", colClasses = "character")
cohort <- melt(cohort, id.vars = c("patient_id", "vih_dt"), value.name = "control_id")
cohort <- cohort[control_id != ""]

# Create object needed after preprocessing:
diagnosis_full <- data.table()
bdu <- data.table()

# Map ICPC codes to ICD9 codes
for (i in c("control/diagnosticos_ctrl", "cohorte/diagnosticos")) {
  
  message(paste0("\nMapping to ICD9: ", i, "\n---------------------------------------------------------------------\n"))
  diagnosis <- suppressWarnings(fread(paste0(i, ".csv"), colClasses = "character", encoding = "UTF-8"))
  diagnosis <-  map_cie9(diagnosis = diagnosis, 
                         dictionary = "mapeo_ciap_cie9/diccionario_ciapAragon_bifap_cie9.csv", 
                         wrong_descriptors = "mapeo_ciap_cie9/1_corregir_descriptores_originales.csv")
  
  fwrite(diagnosis[!is.na(CIE9)], paste0(i, "_cie9.csv"), encoding = "UTF-8")
  
}

# Transformn diagnosis to ccs:
for (i in c("control/", "cohorte/")) {
  diagnosis <- fread(paste0(i,grep("_cie9.csv",list.files(i), value = T)), colClasses = "character", encoding = "UTF-8")
  diagnosis[, cie9_no_dot := gsub("\\.", "", CIE9)]
  diagnosis <- merge(diagnosis, cci[,.(cie9_no_dot, chronic_bool)], by = "cie9_no_dot", all.x = T)
  diagnosis <- merge(diagnosis, ccs, by = "cie9_no_dot", all.x = T)
  diagnosis <- diagnosis[!is.na(ccs_epichron_code)]
  
  setkey(diagnosis,NULL)
  # Set date variables to proper data type:
  diagnosis[, vih_dt := as.Date(vih_dt)]
  diagnosis[, diag_dt := as.Date(diag_dt)]
  # Set order and drop duplicates
  setorder(diagnosis, patient_id, vih_dt, ccs_epichron_code, diag_dt)
  diagnosis <- unique(diagnosis[chronic_bool == 1], by = c("patient_id", "vih_dt","ccs_epichron_code"))
  # diagnosis[, ccs_category := paste0(ccs_epichron_code, "_", ccs_epichron_label)]

  # Save full information
  setkey(t, NULL)
  diagnosis_full <- rbindlist(list(diagnosis_full, diagnosis), fill = T)
  
  rm(diagnosis)
  
  # Create unified bdu
  temp <- fread(paste0(i, grep("bdu",list.files(i), value = T)), encoding = "UTF-8")
  temp[, year := year(nacimiento_dt)]
  temp[, edad := 2022-year]
  bdu <- rbind(bdu, temp)
  rm(temp)
}

# Mark vih+ patient in all datasets
bdu[, vih_bool := patient_id %in% cohort$patient_id]
bdu <- bdu[patient_id %in% cohort[, c(patient_id, control_id)]]

# Save all transform data:
fwrite(bdu, "intermediate/bdu_full.csv", encoding = "UTF-8")
fwrite(diagnosis_full, "intermediate/diagnosis_chronic_full.csv", encoding = "UTF-8")

