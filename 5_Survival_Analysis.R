# Author: Alejandro Santos Mejías
# Date last update: 2025-09-23
# Input: Date incidence diagnosis datasets
# Output: 
# Motivational comment: May death finds you rested.
# 
#==============#
#= Acknoledge =#
#==============#
# Tutorials: 
# https://www.epirhandbook.com/en/new_pages/survival_analysis.html
# https://www.emilyzabor.com/survival-analysis-in-r.html


rm(list = ls())
gc()

source("99_paths_and_packages.R")

# =========================== #
# === Survival analysis ===== #
# =========================== #

# ============================================================================ #

bdu <- fread("intermediate/bdu_full.csv", encoding = "UTF-8")
survival <- fread("intermediate/diagnosis_date_incident_final.csv", encoding = "UTF-8")
cohort <- fread("intermediate/cohort_match_id.csv", encoding = "UTF-8")

# Transformn to date
survival[, HIV_infection := NULL]
cols <- colnames(survival)[-c(1,2)]
survival[, (cols) := lapply(.SD, lubridate::ymd), .SDcols = cols]
survival[, vih_dt := lubridate::ymd(vih_dt)]

# Assign match id and duplicate controls with different vih_dt
cohort[, variable := NULL]
cohort <- melt(cohort, id.vars = c("vih_dt", "match_id"), value.name = "patient_id")
cohort <- unique(cohort, by = c("patient_id", "match_id"))
# Duplicate controls in diagnosis dataset
# cohort[, .N, by = "patient_id"] # Check which one
survivalClean <- merge(survival, cohort[, .(patient_id, match_id)], by = "patient_id")
setcolorder(survivalClean, "match_id")
rm(survival)
# Update vih_dt
survivalClean[cohort[, .(match_id, vih_dt)], on = .(match_id), vih_dt := i.vih_dt]

# Prepare dataset for survival analysis
survivalAnal <- melt(survivalClean, id.vars = c("match_id", "patient_id", "vih_dt"), variable.name = "disease", value.name = "disease_dt")
# Check no disease_dt before starting point vih_dt
survivalAnal[vih_dt > disease_dt, disease_dt := NA]
idNopDis <- survivalAnal[is.na(disease_dt), .N, by = "patient_id"][N == length(cols), patient_id]
# Extract all censored patient_id + vih_dt registries with no diseases previous end of study 2022-12-31
idNopDis <- unique(survivalAnal[patient_id %in% idNopDis, .(match_id, patient_id, vih_dt)])
idNopDis <- merge(idNopDis, bdu[, .(patient_id, sexo, edad, baja_dt, vih_bool)], by = "patient_id")
idNopDis[, baja_dt := as.Date(baja_dt)]
idNopDis[vih_dt >= baja_dt, baja_dt := vih_dt]
idNopDis[is.na(baja_dt), baja_dt := as.Date("20221231", format = "%Y%m%d")]
setnames(idNopDis, "baja_dt", "disease_dt")

# Extract all non-censored patient_id + vih_dt 
survivalAnal <- survivalAnal[!is.na(disease_dt)]
survivalAnal[, patient_idvih_dt := paste0(patient_id, vih_dt)]
# Extract by index those registers with the first disease_dt
ind <- survivalAnal[, .(row_index = .I[which.min(disease_dt)]), by = "patient_idvih_dt"][, row_index]
survivalAnal <- survivalAnal[ind]
survivalAnal[, patient_idvih_dt := NULL]
# Create final dataset
survivalAnal <- merge(survivalAnal, bdu[, .(patient_id, sexo, edad, vih_bool)], by = "patient_id")
survivalAnal <- rbindlist(list(survivalAnal, idNopDis), fill = T)
# Censor status and times in days to the event 
# censor status 0 = censored, 1 = event happening, remember to code it as numeric, survival is extremely picky with data type
survivalAnal[, `:=`(time = as.numeric(disease_dt - vih_dt)/365.25,
                    censored_status = ifelse(is.na(disease), 0, 1))]
# Prepared rest of groups:
age_cuts <- c(-Inf,44,65, Inf)
age_labels <- c("under 45", "45 to 65", "over 65")
survivalAnal[, agebands := cut(edad, age_cuts, labels = age_labels)]
survivalAnal[,  `:=`(vih_bool = factor(vih_bool, levels = c(F,T), labels = c("NON-HIV", "HIV")),
                     sexo = factor(sexo, levels = c("HOMBRE", "MUJER"), labels = c("Men", "Women")))]
# Filter out patient censored at time 0
survivalAnal <- survivalAnal[ !(time == 0 & censored_status == 1)]
# Run the model
survSex <- survfit2(Surv(time = time, event = censored_status) ~ vih_bool, data = survivalAnal)

# Plot the results stratified by ageband and sex
p <- survminer::ggsurvplot(
  survSex, 
  data = survivalAnal,          # again specify the data used to fit linelistsurv_fit_sex 
  conf.int = T,              # show confidence interval of KM estimates
  surv.scale = "percent",        # present probabilities in the y axis in %
  break.time.by = 5,            # present the time axis with an increment of 10 days
  xlab = "Follow-up years",
  ylab = "Survival Probability",
  facet.by = c("agebands", "sexo"),
  censor.shape = "",
  pval = T,                      # print p-value of Log-rank test
  # pval.coord = c(100,10),        # print p-value at these plot coordinates
  # risk.table = T,                # print the risk table at bottom 
  legend.title = "HIV status: ",       # legend characteristics
  # legend.labs = c("Female","Male"),
  # font.legend = 10, 
  palette = c("#f58231", "#7E0CD6", "#1cb7fa"),             # specify color palette
  # palette = "Dark",
  surv.median.line = "hv",       # draw horizontal and vertical lines to the median survivals
  ggtheme = theme_light(),        # simplify plot background
  panel.labs = list(sexo = c("Men", "Women"), agebands = c("Under 45 years", "Between 45 to 65 years", "Over 65 years")),
  short.panel.labs = T
)

p + theme(strip.text = element_text(size = 12, face = "bold", colour = "black"), legend.title = element_text(face = "bold"))

ggsave("paper/survival/Survival_curves.png", width = 2250, height = 2625, unit = "px", dpi = 300)
ggsave("paper/survival/Survival_curves.tiff", width = 2250, height = 2625, unit = "px", dpi = 300)

survSex <- survfit2(Surv(time = time, event = censored_status) ~ vih_bool + sexo + agebands, data = survivalAnal)
write(capture.output(survSex), "paper/survival/medianSurvivalTime.txt")

#### Top 10 most adquired pathologies: ####
t <- survivalAnal[, .(n_grupo = .N), by = c("vih_bool", "sexo", "agebands", "disease")][order(n_grupo)]
t <- t[survivalAnal[, .(Total_grupo = .N), by = c("vih_bool", "sexo", "agebands")], on = c("vih_bool", "sexo", "agebands")]
t[, proporcion := round(n_grupo/Total_grupo * 100, 3)]
fwrite(t, file = "paper/survival/descriptivoSurv.csv", encoding = "UTF-8")
setorder(t, -proporcion)
n <- 10
fwrite(t[, .(diseses = head(disease, n), n_group = head(n_grupo, n), Total_grupo = head(Total_grupo,n) ,prop_disease = head(proporcion, n)), by = c("vih_bool", "sexo", "agebands")], file = "paper/survival/desciptivoTopN.csv")


#### Sensitivity analysis with homogeneous number of diseases at the moment of matching ####
t <- fread("intermediate/diagnosis_date_prevalent_final.csv", encoding = "UTF-8", colClasses = "character")
t[, HIV_infection := NULL]
cols <- colnames(t)[-c(1,2)]
t[, (cols) := lapply(.SD, function(x){as.Date(as.character(x), format = "%Y%m%d")}), .SDcols = cols]
t[, vih_dt := as.Date(vih_dt)]
t[,(cols) := lapply(.SD, function(x){ifelse(vih_dt < x | is.na(x), 0, 1)}), .SDcols = cols]
t[, nDis := rowSums(.SD), .SDcols = cols]
cohort <- merge(cohort, t[, .(patient_id, nDis)], by = "patient_id")
cohort[, diff_nDis := abs(nDis - nDis[variable == "patient_id"]), by = match_id]

# Distribution:
cohort[, .N, by = diff_nDis][order(diff_nDis)]

t <- survivalAnal[patient_id %in% cohort[diff_nDis <= 0, patient_id]]
# Only people with a pair:
t <- t[match_id %in% t[, .N, by = match_id][N >= 3, match_id]]
# t <- unique(t, by = c("match_id", "vih_bool"))

survSex <- survfit2(Surv(time = time, event = censored_status) ~ vih_bool, data = t)

p <- survminer::ggsurvplot(
  survSex, 
  data = t,          # again specify the data used to fit linelistsurv_fit_sex 
  conf.int = T,              # show confidence interval of KM estimates
  surv.scale = "percent",        # present probabilities in the y axis in %
  break.time.by = 5,            # present the time axis with an increment of 10 days
  xlab = "Follow-up years",
  ylab = "Survival Probability",
  facet.by = c("agebands", "sexo"),
  censor.shape = "",
  pval = T,                      # print p-value of Log-rank test
  # pval.coord = c(100,10),        # print p-value at these plot coordinates
  # risk.table = T,                # print the risk table at bottom 
  legend.title = "HIV status: ",       # legend characteristics
  # legend.labs = c("Female","Male"),
  # font.legend = 10, 
  palette = c("#f58231", "#7E0CD6", "#1cb7fa"),             # specify color palette
  # palette = "Dark",
  surv.median.line = "hv",       # draw horizontal and vertical lines to the median survivals
  ggtheme = theme_light(),        # simplify plot background
  panel.labs = list(sexo = c("Men", "Women"), agebands = c("Under 45 years", "Between 45 to 65 years", "Over 65 years")),
  short.panel.labs = T
)

p + theme(strip.text = element_text(size = 12, face = "bold", colour = "black"), legend.title = element_text(face = "bold"))

ggsave("paper/survival/Survival_curvesAdjustedByNDis.png", width = 2250, height = 2625, unit = "px", dpi = 300)
ggsave("paper/survival/Survival_curvesAdjustedByNDis.tiff", width = 2250, height = 2625, unit = "px", dpi = 300)


survSex <- survfit2(Surv(time = time, event = censored_status) ~ vih_bool + sexo + agebands, data = t)
write(capture.output(survSex), "paper/survival/medianSurvivalTimeAdjustedByNDis.txt")

#### Top 10 most adquired pathologies in homogeneous multimorbidity burden cohort: ####
t <- survivalAnal[, .(n_grupo = .N), by = c("vih_bool", "sexo", "agebands", "disease")][order(n_grupo)]
t <- t[survivalAnal[, .(Total_grupo = .N), by = c("vih_bool", "sexo", "agebands")], on = c("vih_bool", "sexo", "agebands")]
t[, proporcion := round(n_grupo/Total_grupo * 100, 3)]
fwrite(t, file = "paper/survival/descriptivoSurvAdjustedByNDis.csv", encoding = "UTF-8")
setorder(t, -proporcion)
n <- 10
fwrite(t[, .(diseses = head(disease, n), n_group = head(n_grupo, n), Total_grupo = head(Total_grupo,n) ,prop_disease = head(proporcion, n)), by = c("vih_bool", "sexo", "agebands")], file = "paper/survival/desciptivoTopNAdjustedByNDis.csv")
