library(gtsummary)
library(data.table)
library(lubridate)
library(stringr)
library(ggplot2)
library(openxlsx)


setwd("C:/Santos/VIH/Prueba_push/data7/")
options(scipen = 999)

source("mapeo_ciap_cie9/mapping_ciap_bifap_cie9_v3.0.R")


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
diag_bool_prev <- data.table()
diag_date_prev <- data.table()
diag_bool_inc <- data.table()
diag_date_inc <- data.table()
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
  
  # Get prevalent cases in wide boolean format:
  t_bool_prev <- data.table::dcast(diagnosis[, .(patient_id, vih_dt, ccs_vih_label, chronic_bool)], formula = patient_id + vih_dt ~ ccs_vih_label, value.var = "chronic_bool")
  # Get incident cases in wide boolena format:
  t_bool_inc <- data.table::dcast(diagnosis[vih_dt < diag_dt, .(patient_id, vih_dt, ccs_vih_label, chronic_bool)], formula = patient_id + vih_dt ~ ccs_vih_label, value.var = "chronic_bool")
  
  # Repeat the process for date format:
  diagnosis[, diag_dt := format(diag_dt, "%Y%m%d")]
  t_date_prev <- data.table::dcast(diagnosis[, .(patient_id, vih_dt, ccs_vih_label, diag_dt)], formula = patient_id + vih_dt ~ ccs_vih_label, value.var = "diag_dt")
  t_date_inc <- data.table::dcast(diagnosis[ymd(vih_dt) < ymd(diag_dt), .(patient_id, vih_dt, ccs_vih_label, diag_dt)], formula = patient_id + vih_dt ~ ccs_vih_label, value.var = "diag_dt")
  
    # Save full information of all transformn diagnosis, diagnosis in boolean form of prevalent and incident chronic diseases:
  setkey(t, NULL)
  diagnosis_full <- rbindlist(list(diagnosis_full, diagnosis), fill = T)
  diag_bool_prev <- rbindlist(list(diag_bool_prev, t_bool_prev), fill = T)
  diag_bool_inc <- rbindlist(list(diag_bool_inc, t_bool_inc), fill = T)
  diag_date_prev <- rbindlist(list(diag_date_prev, t_date_prev), fill = T)
  diag_date_inc <- rbindlist(list(diag_date_inc, t_date_inc), fill = T)
  
  rm(t_bool_prev, t_bool_inc, t_date_prev, t_date_inc, diagnosis)
  
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
diag_bool_prev[patient_id %in% cohort$patient_id, `HIV infection` := 1]
diag_bool_inc[patient_id %in% cohort$patient_id, `HIV infection` := 1]
diag_date_prev[patient_id %in% cohort$patient_id, `HIV infection` := format(vih_dt, "%Y%m%d")]
diag_date_inc[patient_id %in% cohort$patient_id, `HIV infection` := format(vih_dt, "%Y%m%d")]

# Keep only adult subjects:
diag_bool_prev <- diag_bool_prev[patient_id %in% bdu[edad >= 18, patient_id]]
diag_bool_inc <- diag_bool_inc[patient_id %in% bdu[edad >= 18, patient_id]]
diag_date_prev <- diag_date_prev[patient_id %in% bdu[edad >= 18, patient_id]]
diag_date_inc <- diag_date_inc[patient_id %in% bdu[edad >= 18, patient_id]]

# Patients with only one match control, assume they are completely healthy and add those control to the diag_bool datasets:
fixUnevenMatching <- function(data = NA, cohort = NA){
  
  cohort_missing <- cohort[patient_id %in% data[, c(patient_id)]][control_id %in% data[, c(patient_id)]]
  patient_id_1_match <- cohort_missing[, .N, by = patient_id][N == 1, patient_id]
  cohort_missing <- cohort[patient_id %in% patient_id_1_match & !(control_id %in% cohort_missing[, c(patient_id)]), .(control_id, vih_dt)]
  setnames(cohort_missing, c("patient_id", "vih_dt"))
  cohort_missing[, vih_dt := as.Date(vih_dt)]
  data <- rbindlist(list(data, cohort_missing), fill = T)
  return(data)
}

diag_bool_prev <- fixUnevenMatching(data = diag_bool_prev, cohort = cohort)
diag_bool_inc <- fixUnevenMatching(data = diag_bool_inc, cohort = cohort)
diag_date_prev <- fixUnevenMatching(data = diag_date_prev, cohort = cohort)
diag_date_inc <- fixUnevenMatching(data = diag_date_inc, cohort = cohort)

# Set all NA to 0 in boolean datasets:
cols <- colnames(diag_bool_prev)
diag_bool_prev[, (cols[-c(1:2)]) := lapply(.SD, function(x) ifelse(is.na(x),0,1)), .SDcols = cols[-c(1:2)]]

cols <- colnames(diag_bool_inc)
diag_bool_inc[, (cols[-c(1:2)]) := lapply(.SD, function(x) ifelse(is.na(x),0,1)), .SDcols = cols[-c(1:2)]]

# Heal colnames to avoid issues in formulas and reorder:
fixColnamesFormulas <- function(data){
  setnames(data, gsub(",|;|:| |`|\\(|\\)", "_", colnames(data)))
  setcolorder(data, c("patient_id", "vih_dt", "HIV_infection"))
}

fixColnamesFormulas(diag_bool_prev)
fixColnamesFormulas(diag_bool_inc)
fixColnamesFormulas(diag_date_prev)
fixColnamesFormulas(diag_date_inc)


# Save all transform data:
fwrite(bdu, "intermediate/bdu_full.csv", encoding = "UTF-8")
fwrite(diagnosis_full, "intermediate/diagnosis_chronic_full.csv", encoding = "UTF-8")
fwrite(diag_bool_prev,"intermediate/diagnosis_bool_prevalent.csv", encoding = "UTF-8")
fwrite(diag_bool_inc, "intermediate/diagnosis_bool_incident.csv", encoding = "UTF-8")
fwrite(diag_date_prev, "intermediate/diagnosis_date_prevalent.csv", encoding = "UTF-8")
fwrite(diag_date_inc, "intermediate/diagnosis_date_incident.csv", encoding = "UTF-8")
