library(data.table)
library(stringr)
library(ggplot2)

setwd("C:/Santos/VIH/datos7/")

ccs <- fread("ccs_r/ccs_icd9_dic.csv")
# diagnosis <- fread("cohorte/diagnosticos_cie9.csv", colClasses = "character", encoding = "UTF-8")
diagnosis <- fread("control/diagnosticos_ctl_cie9.csv", colClasses = "character", encoding = "UTF-8")

diagnosis[, cie9_no_dot := gsub("\\.", "", CIE9)]

cci <- fread("ccs_r/cci2015.csv", quote = "")
# cci[, `ICD-9-CM CODE` := gsub(" ","",`ICD-9-CM CODE`)]
# cci[, `ICD-9-CM CODE` := gsub("'","",`ICD-9-CM CODE`)]
# cci[, `ICD-9-CM CODE DESCRIPTION` := gsub("\\\"", "", `ICD-9-CM CODE DESCRIPTION`)]
# cci[, `CATEGORY DESCRIPTION` := gsub("'", "", `CATEGORY DESCRIPTION`)]
# cci[, `BODY SYSTEM` := gsub("'", "", `BODY SYSTEM`)]
# setnames(cci, c("cie9_no_dot", "cie9_label", "chronic_bool", "body_system"))
# fwrite(cci, "ccs_r/cci2015.csv")
diagnosis <- merge(diagnosis, cci[,.(cie9_no_dot, chronic_bool)], by = "cie9_no_dot", all.x = T)
# diagnosis <- diagnosis[chronic_bool == 1]
diagnosis <- merge(diagnosis, ccs, by = "cie9_no_dot", all.x = T)

# bdu <- fread("cohorte/bdu.csv", encoding = "UTF-8")
bdu <- fread("control/bdu_ctrl.csv", encoding = "UTF-8")

bdu[, year := year(nacimiento_dt)]
bdu[, edad := 2022-year]

setkey(diagnosis,NULL)
setorder(diagnosis, patient_id, vih_dt, ccs_code, diag_dt)
# unique(diagnosis, by = c("patient_id", "vih_dt", "ccs_cd"))
# summary(merge(diagnosis[CIE9 != "042", .N, by = patient_id], bdu[, .(patient_id, edad)]))
# summary(merge(unique(diagnosis, by = c("patient_id", "vih_dt", "ccs_cd"))[CIE9 != "042", .N, by = patient_id], bdu[, .(patient_id, edad)]))
# summary(merge(unique(diagnosis, by = c("patient_id", "vih_dt", "ccs_cd"))[CIE9 != "042", .N, by = patient_id][N > 1], bdu[, .(patient_id, edad)]))
# 
# summary(merge(diagnosis[chronic_bool == 1 & CIE9 != "042", .N, by = patient_id], bdu[, .(patient_id, edad)]))
# summary(merge(unique(diagnosis, by = c("patient_id", "vih_dt", "ccs_cd"))[chronic_bool == 1 & CIE9 != "042", .N, by = patient_id], bdu[, .(patient_id, edad)]))
# summary(merge(unique(diagnosis, by = c("patient_id", "vih_dt", "ccs_cd"))[chronic_bool == 1 & CIE9 != "042", .N, by = patient_id][N > 1], bdu[, .(patient_id, edad)]))

# Check graphically number of person by number of diagnosis
# t <- as.data.table(table(unique(diagnosis, by = c("patient_id", "vih_dt", "ccs_cd"))[chronic_bool == 1 & CIE9 != "042", .N, by = patient_id][,c(N)]))
# t[, V1 := factor(as.integer(V1))]
# 
# ggplot(t[1:11], (aes(x = V1, y = N))) + geom_bar(stat = "identity") +
#   geom_text(aes(label = N, vjust = -0.5))

# Transform into a binary matrix using mininmal date as reference, only taken chronic conditions for this aim:
diagnosis <- unique(diagnosis[chronic_bool == 1], by = c("patient_id", "vih_dt","ccs_code"))
diagnosis[, ccs_category := paste0(ccs_code, "_", ccs_label)]
t <- reshape(diagnosis[, .(patient_id, diag_dt, ccs_category)], idvar = "patient_id", timevar = "ccs_category", direction = "wide")
colnames(t) <- gsub("diag_dt\\.", "", colnames(t))
cols <- colnames(t)
t[, (cols[-1]) := lapply(.SD, function(x) ifelse(is.na(x),0,1)), .SDcols = cols[-1]]
# t[, `5_HIV infection` := 1]
t[, `5_HIV infection` := 0]
# fwrite(diagnosis, "cohorte/diagnosticos_ccs.csv")
# fwrite(diagnosis, "control/diagnosticos_ctl_ccs.csv")


# t2 <- as.data.table(table(diagnosis$ccs_code))
# temp <- merge(t, t2, by = "V1", all = T)
# setnames(temp, c("ccs_code", "N_ctl", "N_vih"))
# ccs[, ccs_code := as.character(ccs_code)]
# temp <- merge(temp, unique(ccs[,.(ccs_code, ccs_label)]), by = "ccs_code")
