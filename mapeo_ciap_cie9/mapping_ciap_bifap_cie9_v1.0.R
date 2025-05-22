library(data.table)
library(stringr)

setwd("C:/Santos/VIH/datos7/")
n1 <- Sys.time()
diagnosis <- fread("control/diagnosticos_ctrl.csv", encoding = "UTF-8", colClasses = "character")
# Please set colnames of your dataset to the proper one
# patient_id, diag_st, diag_cd will be used for the following tasks of mapping
# TODO a function to automateate this task

#### Preprocessing ####
# Set strings and code to upper case:
cols <- c("diag_cd","diag_st")
diagnosis[, (cols) := lapply(.SD,toupper), .SDcols = cols]

# Save number of original registries
n_original <- nrow(diagnosis)

# Remove registries that are wrong, empty or correspond to none assist 
diagnosis <- diagnosis[patient_id != "" | !is.na(patient_id)] 
diagnosis <- diagnosis[grepl("\\s+ERROR|^ERROR", diag_st) == F]
diagnosis <- diagnosis[grepl("ERRONEO|EQUIVOCA", diag_st) == F]
diagnosis <- diagnosis[!(diag_cd %in% c("SAL", "XXX", "Z31", "Q01"))]
diagnosis <- diagnosis[grepl("NO\\s+ACUDE", diag_st) == F]

# Remove doubtful diagnosis
diagnosis <- diagnosis[grepl("\\?|\\¿", diag_st) == F]

# Correct typos
diagnosis[,diag_cd := gsub("O", "0", diag_cd)]
diagnosis[grepl("TB", diag_cd), diag_cd := "T81"]
diagnosis[grepl("K87", diag_cd), diag_cd := "K8A7"]

n_original - nrow(diagnosis)

#### 1 Heal mistaken descriptors ####
descriptors <- fread("mapeo_ciap_cie9/1_corregir_descriptores_originales.csv", colClasses = "character")
dict <- descriptors$descriptor_bueno
names(dict) <- paste0("^",descriptors$descriptor_malo, "$")
# str_replace_all can subtitute multiple pattern for their corresponding values
# if they are provided in a named vector --> name of the vector(pattern to search) = value of the vector(replacement)
# it works as a regex, so the pattern must fit into one, take care of special characters
diagnosis[, diag_st := str_replace_all(diag_st, pattern = dict)]
# TODO Preguntar bea si los cambios son a coincidencias exactas por ejemplo CEFALEA NO ESPECIFICADA

#### 2 Merge labels with BIFAP labels ####
bifap <- fread("mapeo_ciap_cie9/CIAP_BIFAP.txt", colClasses = "character")
bifap[, ISPREF := NULL]
setnames(bifap, "CODE", "BIFAP_CODE")
temp <- merge(diagnosis, bifap[, .(BIFAP_CODE, STR)], by.x= "diag_st", by.y = "STR")

# Remove this BIFAP codes by the grace and power of Luis and Jonas (clinicians)
temp <- temp[!(BIFAP_CODE %in% c("P03.10","P03.19","P03.20","P03.3","P03.7","P03.8","P29.11","P29.17","P29.18","P99.2"))]

#### 3 Assign to BIFAP codes their standard ICD9 label ####
cie9_label <- fread("mapeo_ciap_cie9/PREFERENTES.txt", colClasses = "character")
setnames(cie9_label, c("CODE", "STR"), c("BIFAP_CODE", "BIFAP_LABEL"))
temp <- merge(temp, cie9_label[, .(BIFAP_CODE, BIFAP_LABEL, ETIQUETA_CIE9)], by = "BIFAP_CODE")

#### 4 Map standard ICD9 labels to their corresponding ICD9 codes and ICPC2 codes ####
# This mapping belong to Navarra ... (insert citation)
cie9_codes <- fread("mapeo_ciap_cie9/CIAP_CIE9.txt", colClasses = "character")
setnames(cie9_codes, "DESCRIPTOR", "ETIQUETA_CIE9")
temp <- merge(temp, cie9_codes, by = "ETIQUETA_CIE9")
# fwrite(temp[, .(patient_id, CIE9, ETIQUETA_CIE9, fec_apert)], "control/diagnosticos_ctl_cie9.csv")




nrow(temp)/n_original
# setorder(diagnosis, patient_id, diag_dt, vih_dt)
# setorder(temp, patient_id, diag_dt, vih_dt)
# diagnosis[!(diagnosis[, paste0(patient_id, diag_dt, vih_dt)] %in% temp[, paste0(patient_id, diag_dt, vih_dt)])]
n2 <- Sys.time()
n2-n1
# TODO check with Bea processing time in comparison with the sql way

# Unified dictionary
bifap <- bifap[!(BIFAP_CODE %in% c("P03.10","P03.19","P03.20","P03.3","P03.7","P03.8","P29.11","P29.17","P29.18","P99.2")), .(BIFAP_CODE, STR)]
cie9_label[, c("ISPREF", "AUI") := NULL]
t <- merge(bifap, cie9_label, by = "BIFAP_CODE", all = T)
t <- merge(t, cie9_codes, by = "ETIQUETA_CIE9", all = T)
fwrite(t, "mapeo_ciap_cie9/diccionario_ciapAragon_bifap_cie9.csv")
temporal <- merge(diagnosis, t, by.x = "diag_st", by.y = "STR")
temporal[!is.na(CIE9), .N]
temporal[!is.na(CIE9), .N]- temp[,.N]
# It seems to work fine as a unified dictionary