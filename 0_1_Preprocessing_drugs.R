library(gtsummary)
library(data.table)
library(lubridate)
library(stringr)
library(ggplot2)
library(openxlsx)


setwd("C:/Santos/VIH/Prueba_push/data7/")
options(scipen = 999)

##############################
##### Preprocessing Drugs ####
##############################

################################################################################
# Load dictionaries:
drugAtc <- fread("data_farmacos/Tabla_ATC_BDCAP_2023.csv", colClasses = "character", encoding = "UTF-8")
drugAtc <- drugAtc[nivel_atc == 5]

drugDisease <- fread("data_farmacos/dicHivDrugs2ChronicDiseases.csv", encoding = "UTF-8", colClasses = "character")

cohort <- fread("control/controles.csv", encoding = "UTF-8", colClasses = "character")
cohort <- melt(cohort, id.vars = c("patient_id", "vih_dt"), value.name = "control_id")
cohort <- cohort[control_id != ""]

# Create object needed after preprocessing:
diag_bool_prev <- data.table()
diag_date_prev <- data.table()
diag_bool_inc <- data.table()
diag_date_inc <- data.table()

# Load drugs prescription dataset:
proxyDrugs <- suppressWarnings(fread("data_farmacos/phar_presc.csv", encoding = "UTF-8", colClasses = "character"))[, .(patient_id, vih_dt, treat_start_dt, pharactivesubs_cd)]
# Properly format variables:
proxyDrugs[, vih_dt := as.Date(vih_dt)]
proxyDrugs[, treat_start_dt := ymd(treat_start_dt)]
setorder(proxyDrugs, patient_id, vih_dt, -treat_start_dt, pharactivesubs_cd)
proxyDrugs <- unique(proxyDrugs)
setnames(proxyDrugs, c("treat_start_dt", "pharactivesubs_cd"), c("diag_dt", "atc_cd"))

# Map ATC code to ATC label:
proxyDrugs <- merge(proxyDrugs, drugAtc[, .(atc_cd, atc_label)], by = "atc_cd")

# Map active pharmaceutical substance to CCS disease category:
proxyDrugs <- merge(proxyDrugs, drugDisease, by.x = c("atc_cd", "atc_label"), by.y= c("atc_code", "drug_label"))
setorder(proxyDrugs, patient_id, vih_dt, ccs_label, diag_dt)
proxyDrugs <- unique(proxyDrugs, by = c("patient_id", "vih_dt", "ccs_label"))

# Heal CCS names to avoid issues in formulas:
proxyDrugs[, ccs_label := gsub(",|;|:| |`|\\(|\\)", "_", ccs_label)]

# Change diag_dt format to string:
proxyDrugs[, diag_dt:= format(diag_dt, "%Y%m%d")]

# Save data set for future checks:
fwrite(proxyDrugs, file = "intermediate/diagnosis_chronic_drugs_full.csv", encoding = "UTF-8")

# Create unified enriched diagnoses dataset from both data sources:
diag <- fread("intermediate/diagnosis_chronic_full.csv", encoding = "UTF-8", colClasses = "character")
diag <- diag[, .(patient_id, vih_dt, diag_dt, ccs_vih_label, ccs_epichron_code, cie9_no_dot, ETIQUETA_CIE9)]
# Create an source origin variable:
diag[, source_origin := "Diagnoses"]
# Set same colnames:
setnames(diag, 
         c("patient_id" ,"vih_dt", "diag_dt", "ccs_vih_label", "ccs_epichron_code", "cie9_no_dot", "ETIQUETA_CIE9", "source_origin"),
         c("patient_id", "vih_dt", "diag_dt", "ccs_label", "ccs_code", "source_code", "source_label", "source_origin"))
setcolorder(diag,  c("patient_id", "vih_dt", "diag_dt" ,"ccs_label", "ccs_code", "source_code", "source_label", "source_origin"))
diag[, ccs_label := gsub(",|;|:| |`|\\(|\\)", "_", ccs_label)]
diag[, vih_dt := as.Date(vih_dt)]
proxyDrugs[, source_origin := "Drugs"]
setnames(proxyDrugs, 
         c("patient_id", "vih_dt", "diag_dt" ,"ccs_label", "ccs_code", "atc_cd", "atc_label", "source_origin"),
         c("patient_id", "vih_dt", "diag_dt" ,"ccs_label", "ccs_code", "source_code", "source_label", "source_origin"))
setcolorder(proxyDrugs,  c("patient_id", "vih_dt", "diag_dt" ,"ccs_label", "ccs_code", "source_code", "source_label", "source_origin"))
# Bind together and remove newer registry for each patient and ccs category:
ccsFull <- rbindlist(list(proxyDrugs, diag))
ccsFull[, diag_dt := ymd(diag_dt)]
setorder(ccsFull, patient_id, vih_dt, diag_dt, ccs_label, source_origin)
ccsFull <- unique(ccsFull, by = c("patient_id", "vih_dt", "ccs_label"))


# Remove people present in cohort that are missing in BDU:
bdu <- fread("intermediate/bdu_full.csv", colClasses = "character", encoding = "UTF-8")
cohort <- cohort[patient_id %in% bdu$patient_id & control_id %in% bdu$patient_id]
cohort[, vih_dt := as.Date(vih_dt)]
fwrite(cohort,"intermediate/cohort_final.csv", encoding = "UTF-8")
bdu <- bdu[patient_id %in% cohort[, c(patient_id, control_id)]]
bdu <- bdu[edad >= 18]
fwrite(bdu, "intermediate/bdu_full.csv", encoding = "UTF-8")

# Prevalent cases:
# Switch to wide format:
ccsFull[, diag_dt := format(diag_dt, "%Y%m%d")]
ccsFull <- ccsFull[patient_id %in% bdu$patient_id]
fwrite(ccsFull, file = "intermediate/diagnosis_full_final.csv", encoding = "UTF-8")
diag_date_prev <- data.table::dcast(data = ccsFull[, .(patient_id, vih_dt, ccs_label, diag_dt)], formula = patient_id + vih_dt ~ ccs_label, value.var = "diag_dt")

# Patients with only one match control, assume they are completely healthy and add those control to the datasets:
fixUnevenMatching <- function(data = NA, cohort = NA){
  
  cohort_missing <- cohort[patient_id %in% data$patient_id & !(control_id %in% data$patient_id), .(control_id, vih_dt)]
  setnames(cohort_missing, c("patient_id", "vih_dt"))
  cohort_missing[, vih_dt := as.Date(vih_dt)]
  data <- rbindlist(list(data, cohort_missing), fill = T)
  return(data)
}

diag_date_prev <- fixUnevenMatching(data = diag_date_prev, cohort = cohort)
fwrite(diag_date_prev, file = "intermediate/diagnosis_date_prevalent_final.csv", encoding = "UTF-8")
diag_bool_prev <- copy(diag_date_prev)
cols <- colnames(diag_bool_prev)[-c(1,2)]
diag_bool_prev <- diag_bool_prev[, (cols) := lapply(.SD, function(x){ifelse(is.na(x), 0, 1)}), .SDcols = cols]
fwrite(diag_bool_prev, file = "intermediate/diagnosis_bool_prevalent_final.csv", encoding = "UTF-8")

# Incident cases:
diag_date_inc <- data.table::dcast(data = ccsFull[ymd(vih_dt) <= ymd(diag_dt), .(patient_id, vih_dt, ccs_label, diag_dt)], formula = patient_id + vih_dt ~ ccs_label, value.var = "diag_dt")
diag_date_inc <-  rbindlist(list(diag_date_inc, diag_date_prev[!(paste0(patient_id, vih_dt) %in% paste0(diag_date_inc$patient_id, diag_date_inc$vih_dt)), .(patient_id, vih_dt)]), fill = T)
fwrite(diag_date_prev, file = "intermediate/diagnosis_date_incident_final.csv", encoding = "UTF-8")
diag_bool_inc <- copy(diag_date_inc)
cols <- colnames(diag_bool_inc)[-c(1,2)]
diag_bool_inc <- diag_bool_inc[, (cols) := lapply(.SD, function(x){ifelse(is.na(x), 0, 1)}), .SDcols = cols]
fwrite(diag_bool_inc, file = "intermediate/diagnosis_bool_incident_final.csv", encoding = "UTF-8")

# Fix hiv patients without being positive in ccs HIV Infection:
bdu[patient_id %in% cohort$control_id, vih_bool := F]
bdu[patient_id %in% cohort$patient_id, vih_bool := T]
diag_bool_prev[patient_id %in% bdu[vih_bool == T, patient_id], HIV_infection := 1]
diag_bool_prev <- diag_bool_prev[patient_id %in% bdu$patient_id,]

diag_bool_inc[patient_id %in% bdu[vih_bool == T, patient_id], HIV_infection := 1]
diag_bool_inc <- diag_bool_inc[patient_id %in% bdu$patient_id]

diag_date_prev[patient_id %in% bdu[vih_bool == T, patient_id] & is.na(HIV_infection), HIV_infection := format(vih_dt, "%Y%m%d")]
diag_date_prev <- diag_date_prev[patient_id %in% bdu$patient_id]

diag_date_inc[patient_id %in% bdu[vih_bool == T, patient_id] & is.na(HIV_infection), HIV_infection := format(vih_dt, "%Y%m%d")]
diag_date_inc <- diag_date_inc[patient_id %in% bdu$patient_id]

# Set column order and save:
setcolorder(diag_bool_prev, c("patient_id", "vih_dt", "HIV_infection"))
setcolorder(diag_bool_inc, c("patient_id", "vih_dt", "HIV_infection"))
setcolorder(diag_date_prev, c("patient_id", "vih_dt", "HIV_infection"))
setcolorder(diag_date_inc, c("patient_id", "vih_dt", "HIV_infection"))

fwrite(diag_bool_prev, file = "intermediate/diagnosis_bool_prevalent_final.csv", encoding = "UTF-8")
fwrite(diag_bool_inc, file = "intermediate/diagnosis_bool_incident_final.csv", encoding = "UTF-8")
fwrite(diag_date_prev, file = "intermediate/diagnosis_date_incident_final.csv", encoding = "UTF-8")
fwrite(diag_date_inc, file = "intermediate/diagnosis_date_incident_final.csv", encoding = "UTF-8")
