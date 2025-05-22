library(data.table)
library(stringr)

map_cie9 <- function(diagnosis = NA, 
                     diag_cd = "diag_cd", diag_st = "diag_st", # The colname of your original CIAP codes and labels
                     dictionary = "diccionario_ciapAragon_bifap_cie9.csv", # The location of the unified ciap_bifap_cie9 dictionary
                     wrong_descriptors = "1_corregir_descriptores_originales.csv", # The location of the typos dictionary
                     missing_codes = T # In case you do not want none 
                     ){
  
  old <- c(diag_cd, diag_st)
  new <- c("diag_cd", "diag_st")
  setnames(diagnosis, old, new)
  
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
  
  #### 1 Heal mistaken descriptors ####
  descriptors <- fread(wrong_descriptors, colClasses = "character")
  dict <- descriptors$descriptor_bueno
  names(dict) <- paste0("^",descriptors$descriptor_malo, "$")
  diagnosis[, diag_st := str_replace_all(diag_st, pattern = dict)]
  
  #### 2 Map to CIE9 ####
  cie9 <- fread(dictionary, colClasses = "character") 
  temporal <- merge(diagnosis, t, by.x = "diag_st", by.y = "STR")
  
  if (missing_codes == F) {
    temporal <- temporal[!is.na(CIE9)]
  }
  
  return(temporal)
}

setwd("C:/Santos/VIH/datos7/mapeo_ciap_cie9/")

diagnosis <- fread("../control/diagnosticos_ctrl.csv", encoding = "UTF-8", colClasses = "character")
temporal <- map_cie9(diagnosis = diagnosis)
fwrite(temporal[!is.na(CIE9)], "../control/diagnosticos_ctl_cie9.csv")

diagnosis <- fread("../cohorte/diagnosticos.csv", encoding = "UTF-8", colClasses = "character")
temporal <- map_cie9(diagnosis = diagnosis)
fwrite(temporal[!is.na(CIE9)], "../cohorte/diagnosticos_cie9.csv")
