# Author: Alejandro Santos Mejías
# Date last update: 2025-08-28
# Input: ICD9 diagnosis datset
# Output: Frequncy of each ICD9 code in each CCS category
# Motivational comment: Everything passess. 

source("99_paths_and_packages.R", encoding = "UTF-8")

#######################################################
##### Frequency of ICD9 codes in each CCS category ####
#######################################################

################################################################################

# Create output path:
pathOutput <- paste0("output/", format(Sys.Date(),"%Y%m%d"), "/Frequency_Analysis")
if(!dir.exists(pathOutput)){dir.create(pathOutput , recursive = T)}

# Read the whole population
bdu <- fread("intermediate/bdu_full.csv", encoding = "UTF-8")

# Load ccs 
diag_prev <- fread("intermediate/diagnosis_chronic_full.csv")
diag_prev <- diag_prev[patient_id %in% bdu$patient_id,]

diag_prev <- diag_prev[, .(patient_id, vih_dt, diag_cd, diag_st, ccs_epichron_code, ccs_vih_label)]
fwrite(diag_prev[, .N, by= c("ccs_vih_label", "diag_cd", "diag_st")][order(ccs_vih_label, -N)], file = paste0(pathOutput, "/General_frequency.csv"))

# Stratified counting
diag_prev <- merge(diag_prev, bdu[, .(patient_id, sexo, edad)], by = "patient_id")
diag_prev[ , ageband := cut(edad, breaks = c(-Inf, 44, 65, Inf), labels = c("under 45", "45 to 65", "over 65"))]
fwrite(diag_prev[, .N, by = c("sexo", "ageband","ccs_vih_label", "diag_st")][order(sexo, ageband, ccs_vih_label, -N)], file = paste0(pathOutput, "/Stratified_frequency_by_icpc2_label.csv"))
fwrite(diag_prev[, .N, by = c("sexo", "ageband","ccs_vih_label", "diag_cd")][order(sexo, ageband, ccs_vih_label, -N)], file = paste0(pathOutput, "/Stratified_frequency_by_icpc2_code.csv"))
