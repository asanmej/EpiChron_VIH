# Author: Alejandro Santos Mejías
# Date last update: 2025-04-25
# Input: bdu and diagnosis dataset
# Output: Remove from datasets all patients without proper age matching
# Comment: This is an optional script, just test results with proper matching performed.


rm(list = ls())
gc()

source("99_paths_and_packages.R", encoding = "UTF-8")

bdu <- fread("intermediate/bdu_full.csv", encoding = "UTF-8")
diag_bool_prev <- fread("intermediate/diagnosis_bool_prevalent_final.csv", encoding = "UTF-8")
diag_bool_inc <- fread("intermediate/diagnosis_bool_incident_final.csv", encoding = "UTF-8")

cohort <- fread("control/controles.csv", encoding = "UTF-8", colClasses = "character")
cohort <- melt(cohort, id.vars = c("patient_id", "vih_dt"), value.name = "control_id")
cohort <- merge(cohort, bdu[, .(patient_id, year)], by = "patient_id")
setnames(cohort, "year", "year_vih")
cohort <- merge(cohort, bdu[, .(patient_id, year)], by.x = "control_id", by.y = "patient_id")
setnames(cohort, "year", "year_ctl")
cohort[, age_difference := abs(year_ctl - year_vih)]
cohort <- cohort[age_difference < 3]

bdu <- bdu[patient_id %in% unique(cohort[, c(patient_id, control_id)])]
diag_bool_inc <- diag_bool_inc[patient_id %in% bdu$patient_id]
diag_bool_prev <- diag_bool_prev[patient_id %in% bdu$patient_id]
fwrite(bdu, "intermediate/bdu_full.csv", encoding = "UTF-8")
fwrite(diag_bool_prev, file = "intermediate/diagnosis_bool_prevalent_final.csv", encoding = "UTF-8")
fwrite(diag_bool_inc, file = "intermediate/diagnosis_bool_incident_final.csv", encoding = "UTF-8")
fwrite(cohort, "intermediate/cohorte_2000.csv", encoding = "UTF-8")
