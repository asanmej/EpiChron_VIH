# Author: Alejandro Santos Mejías
# Date of last update: 2024-10-31


# -----------------------------------------------------------------------------
library(data.table)
library(stringr)

map_cie9 <- function(diagnosis = NA,
                     patient_id = "patient_id", # The colname of the identifier of the patients
                     diag_cd = "diag_cd", diag_st = "diag_st", # The colname of your original CIAP codes and labels
                     dictionary = "diccionario_ciapAragon_bifap_cie9.csv", # The location of the unified ciap_bifap_cie9 dictionary
                     wrong_descriptors = "1_corregir_descriptores_originales.csv", # The location of the typos dictionary
                     missing_codes = T # In case you do not want none mapped codes 
                     ){
  isdt <- is.data.table(diagnosis)
  
  if(!isdt){
    setDT(diagnosis)
  }
  
  old <- c(patient_id, diag_cd, diag_st)
  new <- c("patient_id","diag_cd", "diag_st")
  setnames(diagnosis, old, new)
  
  #### Preprocessing ####
  # Set strings and code to upper case:
  cols <- c("diag_cd","diag_st")
  diagnosis[, (cols) := lapply(.SD,toupper), .SDcols = cols]
  diagnosis[, (cols) := lapply(.SD,function(x){ifelse(x == "" | grepl("^\\s+$",x), NA, x)}), .SDcols = cols]
  
  
  # Save number of original registries
  n_original <- nrow(diagnosis)
  
  # Remove registries that are wrong, empty or correspond to none assist 
  diagnosis <- diagnosis[patient_id != "" | !is.na(patient_id)] 
  diagnosis <- diagnosis[grepl("\\s+ERROR|^ERROR", diag_st) == F]
  diagnosis <- diagnosis[grepl("ERRONEO|EQUIVOCA", diag_st) == F]
  diagnosis <- diagnosis[!(diag_cd %in% c("SAL", "XXX", "Z31", "Q01"))]
  diagnosis <- diagnosis[grepl("NO\\s+ACUDE", diag_st) == F]
  diagnosis <- diagnosis[!is.na(diag_st)]
  
  # Remove doubtful diagnosis
  diagnosis <- diagnosis[grepl("\\?|\\¿", diag_st) == F]
  
  # Correct typos
  diagnosis[,diag_cd := gsub("O", "0", diag_cd)]
  diagnosis[grepl("TB", diag_cd), diag_cd := "T81"]
  diagnosis[grepl("K8A7", diag_cd), diag_cd := "K87"]
  message("Preprocesamiento terminado")
  message(paste0("Se han eliminado ", n_original-diagnosis[,.N], " registros erroneos"))
  
  #### 1 Heal mistaken descriptors ####
  descriptors <- fread(wrong_descriptors, colClasses = "character", encoding = "UTF-8")
  diagnosis <- merge(diagnosis, descriptors, by.x = "diag_st", by.y = "descriptor_malo", all.x = T)
  diagnosis[!is.na(descriptor_bueno), diag_st := descriptor_bueno]
  .Last.updated
  diagnosis[, descriptor_bueno := NULL]
  message("Etiquetas defectuosas corregidas")
  rm(descriptors)
  gc()
  
  #### 2 Map to CIE9 ####
  cie9 <- fread(dictionary, colClasses = "character", encoding = "UTF-8") 
  temporal <- merge(diagnosis, cie9, by.x = "diag_st", by.y = "STR", all.x = T)
  
  if (missing_codes == F) {
    temporal <- temporal[!is.na(CIE9)]
  }
  
  if(!isdt){
    setDF(temporal)
  }
  message("Mapeo finalizado\n",round(temporal[!is.na(CIE9),.N]/temporal[,.N]*100, 2)," transformadas con exito", "\nPasa buen dia.")
  return(temporal)
}

# setwd("Y:/PROYECTOS/2023 ROC18_SAFETY-VAC/Desarrollo/ETL/mapeo_ciap_cie/")
# 
# 
# diagnosis <- fread("../../Datos2/Raw_data/Diagnosticos_hasta_20231231.csv", encoding = "UTF-8", colClasses = "character")
# t1 <- Sys.time()
# temporal <- map_cie9(diagnosis = diagnosis,
#                      patient_id = "cia_anonimizado",
#                      diag_cd = "codigo_CIAP", 
#                      diag_st = "descriptor")
# t2 <- Sys.time()
# t2-t1
# fwrite(temporal[!is.na(CIE9)], "../../Intermediate/diagnosticos_cie9.csv")
# fwrite(temporal[is.na(CIE9)], "../../Intermediate/diagnosticos_ciap.csv")

