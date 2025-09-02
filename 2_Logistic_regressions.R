# Author: Alejandro Santos Mejías
# Date last update: 2025-04-25
# Input: Boolean diagnosis datasets
# Output: Logistic regression output
# Motivational comment: We are one in a multitude, multitude should be one with us. 

rm(list = ls())
gc()

source("99_paths_and_packages.R")

###############################
##### Logistic regressions ####
###############################

################################################################################

bdu <- fread("intermediate/bdu_full.csv")
diag_bool_prev <- fread("intermediate/diagnosis_bool_prevalent_final.csv")
diag_bool_inc <- fread("intermediate/diagnosis_bool_incident_final.csv")

diag_bool_prev <- merge(diag_bool_prev, bdu[, .(patient_id, edad, sexo)], by = "patient_id")
diag_bool_inc <- merge(diag_bool_inc, bdu[, .(patient_id, edad, sexo)], by = "patient_id")
setcolorder(diag_bool_prev, c("patient_id", "vih_dt", "HIV_infection", "sexo", "edad"))
setcolorder(diag_bool_inc, c("patient_id", "vih_dt", "HIV_infection", "sexo", "edad"))
setkey(diag_bool_inc, NULL)
setkey(diag_bool_prev, NULL)

if(!dir.exists(paste0("output/", format(Sys.Date(),"%Y%m%d")))){dir.create(paste0("output/", format(Sys.Date(),"%Y%m%d")))}

##### Modelling prevalent and incident diseases adjusted by HIV status ####
for (l in c("diag_bool_prev", "diag_bool_inc")) {
  results <- data.table()
  for (i in 6:length(colnames(get(l)))) {
    tryCatch(t <- glm(formula = as.formula(paste0("`",colnames(get(l))[i],"`","~`HIV_infection` + edad")), family = "binomial", data = get(l)),
             warning = function(w){message(paste("no convergence at", l, colnames(get(l))[i]))}
    )
    cols <- c(colnames(get(l))[i], "HIV_infection")
    r <- data.table(ccs_epichron_label = colnames(get(l))[i],
                    n_positive_total = table(get(l)[, ..i])[2],
                    n_positive_ctl = ifelse(length(get(l)[, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 0, N]) == 0, 
                                            yes = 0, 
                                            no = get(l)[, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 0, N]),
                    n_positive_hiv = ifelse(length(get(l)[, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 1, N]) == 0, 
                                            yes = 0, 
                                            no = get(l)[, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 1, N]),
                    n_negative_total = table(get(l)[, ..i])[1],
                    n_negative_ctl = ifelse(length(get(l)[, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 0, N]) == 0, 
                                            yes = 0, 
                                            no = get(l)[, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 0, N]),
                    n_negative_hiv = ifelse(length(get(l)[, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 1, N]) == 0, 
                                            yes = 0, 
                                            no = get(l)[, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 1, N]),
                    # Odd ratio of the variable without intercept
                    OR_HIV = exp(coef(t))["HIV_infection"],
                    # Confidence interval by Wald method
                    lb_HIV = ifelse(as.data.frame(exp(confint.default(t)))["2.5 %"]["HIV_infection",] == Inf, NA, as.data.frame(exp(confint.default(t)))["2.5 %"]["HIV_infection",]),
                    ub_HIV = ifelse(as.data.frame(exp(confint.default(t)))["97.5 %"]["HIV_infection",] == Inf, NA, as.data.frame(exp(confint.default(t)))["97.5 %"]["HIV_infection",]),
                    p_value_HIV = as.data.frame(summary(t)$coefficient)["Pr(>|z|)"]["HIV_infection",],
                    # Odd ratio of edad
                    OR_edad = exp(coef(t))["edad"],
                    # Confidence interval by Wald method
                    lb_edad = ifelse(as.data.frame(exp(confint.default(t)))["2.5 %"]["edad",] == Inf, NA, as.data.frame(exp(confint.default(t)))["2.5 %"]["edad",]),
                    ub_edad = ifelse(as.data.frame(exp(confint.default(t)))["97.5 %"]["edad",] == Inf, NA, as.data.frame(exp(confint.default(t)))["97.5 %"]["edad",]),
                    p_value_edad = as.data.frame(summary(t)$coefficient)["Pr(>|z|)"]["edad",]
    )
    r[is.na(ub_HIV) | is.na(lb_HIV), OR_HIV := NA]
    r[is.na(ub_edad) | is.na(lb_edad), OR_edad := NA]
    results <- rbind(results, r)
  }
  # Remove empty row, round results and mark p values < 0.05:
  results <- results[!is.na(ccs_epichron_label)]
  
  cols <- grep("^n_positive_total|^OR|^lb|^ub|^p_value", colnames(results), value = T)
  results[is.na(n_positive_total)|n_positive_total == 0, (cols) := NA]
  
  cols <- grep("^OR|^lb|^ub|^p_value", colnames(results), value = T)
  results[, (cols) := lapply(.SD, function(x) round(x, digits = 3)), .SDcols = cols]
  
  results[, p_05_HIV := p_value_HIV < 0.05]
  results[, p_05_edad := p_value_edad < 0.05]
  
  results <- results[!(is.na(n_positive_total) | is.na(OR_HIV))]
  
  fwrite(results, paste0("output/",format(Sys.Date(), "%Y%m%d"),"/", l, "_", format(Sys.Date(), "%Y%m%d"), ".csv"))
}

##### Modelling stratified by sex: ####
bdu <- fread("intermediate/bdu_full.csv")
diag_bool_prev <- fread("intermediate/diagnosis_bool_prevalent_final.csv")
diag_bool_inc <- fread("intermediate/diagnosis_bool_incident_final.csv")

diag_bool_prev <- merge(diag_bool_prev, bdu[, .(patient_id, edad, sexo)], by = "patient_id")
diag_bool_inc <- merge(diag_bool_inc, bdu[, .(patient_id, edad, sexo)], by = "patient_id")
setcolorder(diag_bool_prev, c("patient_id", "vih_dt", "HIV_infection", "sexo", "edad"))
setcolorder(diag_bool_inc, c("patient_id", "vih_dt", "HIV_infection", "sexo", "edad"))
setkey(diag_bool_inc, NULL)
setkey(diag_bool_prev, NULL)


for (s in c("HOMBRE", "MUJER")){
  for (l in c("diag_bool_prev", "diag_bool_inc")) {
    results <- data.table()
    for (i in 6:length(colnames(get(l)))) {
      tryCatch(t <- glm(formula = as.formula(paste0("`",colnames(get(l))[i],"`","~`HIV_infection` + edad")), family = "binomial", data = get(l)[sexo == s]),
               warning = function(w){message(paste0("no convergence at ", l," ", s, " ", colnames(get(l))[i]))}
      )
      cols <- c(colnames(get(l))[i], "HIV_infection")
      r <- data.table(ccs_epichron_label = colnames(get(l))[i],
                      n_positive_total = table(get(l)[sexo == s, ..i])[2],
                      n_positive_ctl = ifelse(length(get(l)[sexo == s, ..cols][, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 0, N]) == 0, 
                                              yes = 0, 
                                              no = get(l)[sexo == s, ..cols][, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 0, N]),
                      n_positive_hiv = ifelse(length(get(l)[sexo == s, ..cols][, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 1, N]) == 0, 
                                              yes = 0, 
                                              no = get(l)[sexo == s, ..cols][, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 1, N]),
                      n_negative_total = table(get(l)[sexo == s, ..i])[1],
                      n_negative_ctl = ifelse(length(get(l)[sexo == s , ..cols][, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 0, N]) == 0, 
                                              yes = 0, 
                                              no = get(l)[sexo == s , ..cols][, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 0, N]),
                      n_negative_hiv = ifelse(length(get(l)[sexo == s , ..cols][, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 1, N]) == 0, 
                                              yes = 0, 
                                              no = get(l)[sexo == s , ..cols][, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 1, N]),
                      # Odd ratio of the variable without intercept
                      OR_HIV = exp(coef(t))["HIV_infection"],
                      # Confidence interval by Wald method
                      lb_HIV = ifelse(as.data.frame(exp(confint.default(t)))["2.5 %"]["HIV_infection",] == Inf, NA, as.data.frame(exp(confint.default(t)))["2.5 %"]["HIV_infection",]),
                      ub_HIV = ifelse(as.data.frame(exp(confint.default(t)))["97.5 %"]["HIV_infection",] == Inf, NA, as.data.frame(exp(confint.default(t)))["97.5 %"]["HIV_infection",]),
                      p_value_HIV = as.data.frame(summary(t)$coefficient)["Pr(>|z|)"]["HIV_infection",],
                      # Odd ratio of edad
                      OR_edad = exp(coef(t))["edad"],
                      # Confidence interval by Wald method
                      lb_edad = ifelse(as.data.frame(exp(confint.default(t)))["2.5 %"]["edad",] == Inf, NA, as.data.frame(exp(confint.default(t)))["2.5 %"]["edad",]),
                      ub_edad = ifelse(as.data.frame(exp(confint.default(t)))["97.5 %"]["edad",] == Inf, NA, as.data.frame(exp(confint.default(t)))["97.5 %"]["edad",]),
                      p_value_edad = as.data.frame(summary(t)$coefficient)["Pr(>|z|)"]["edad",]
      )
      r[is.na(ub_HIV) | is.na(lb_HIV), OR_HIV := NA]
      r[is.na(ub_edad) | is.na(lb_edad), OR_edad := NA]
      results <- rbind(results, r)
    }
    # Remove empty row, round results and mark p values < 0.05:
    results <- results[!is.na(ccs_epichron_label)]
    
    cols <- grep("^n_positive_total|^OR|^lb|^ub|^p_value", colnames(results), value = T)
    results[is.na(n_positive_total) | n_positive_total == 0, (cols) := NA]
    
    cols <- grep("^OR|^lb|^ub|^p_value", colnames(results), value = T)
    results[, (cols) := lapply(.SD, function(x) round(x, digits = 3)), .SDcols = cols]
    
    results[, p_05_HIV := p_value_HIV < 0.05]
    results[, p_05_edad := p_value_edad < 0.05]
    
    results <- results[!(is.na(n_positive_total) | is.na(OR_HIV))]
    
    fwrite(results, paste0("output/",format(Sys.Date(), "%Y%m%d"),"/", l, "_", s, "_", format(Sys.Date(), "%Y%m%d"), ".csv"))
  }
}


##### Modelling stratified by sex and agebands: ####
bdu <- fread("intermediate/bdu_full.csv", encoding = "UTF-8")
diag_bool_prev <- fread("intermediate/diagnosis_bool_prevalent_final.csv", encoding = "UTF-8")
# diag_bool_inc <- fread("intermediate/diagnosis_bool_incident_final.csv")

diag_bool_prev <- merge(diag_bool_prev, bdu[, .(patient_id, edad, sexo, tsi, nacionalidad)], by = "patient_id")
# diag_bool_inc <- merge(diag_bool_inc, bdu[, .(patient_id, edad, sexo)], by = "patient_id")

age_cuts <- c(-Inf,44,65, Inf)
age_labels <- c("<45", "45-65", ">65")
diag_bool_prev[, agebands := cut(edad, age_cuts, labels = age_labels)]
# diag_bool_inc[, agebands := cut(edad, age_cuts, labels = age_labels)]
diag_bool_prev[, edad := NULL]
# diag_bool_inc[, edad := NULL]
setcolorder(diag_bool_prev, c("patient_id", "vih_dt", "HIV_infection", "sexo", "agebands", "tsi", "nacionalidad"))
diag_bool_prev[, nacionalidad := ifelse(nacionalidad == "ESPAÑA", "España", "Extranjero")]
diag_bool_prev[, nacionalidad := as.factor(nacionalidad)]
diag_bool_prev[, HIV_infection := as.factor(HIV_infection)]
diag_bool_prev[, tsi := factor(tsi, levels = c("< 18000", "entre 18000 y 100000", "> 100000", "Farmacia gratuita", "No asegurados", "Mutualistas"))]
# setcolorder(diag_bool_inc, c("patient_id", "vih_dt", "HIV_infection", "sexo", "agebands"))
# setkey(diag_bool_inc, NULL)
setkey(diag_bool_prev, NULL)


for (s in c("HOMBRE", "MUJER")){
  for (l in c("diag_bool_prev"
              # , "diag_bool_inc"
              )) {
    for (age_label in age_labels) {
      results <- data.table()
      for (i in 8:length(colnames(get(l)))) {
        tryCatch(expr = {
          
          t <- glm(formula = as.formula(paste0("`",colnames(get(l))[i],"`","~`HIV_infection` + nacionalidad + tsi")), family = "binomial", data = get(l)[sexo == s & agebands == age_label])
          
          cols <- c(colnames(get(l))[i], "HIV_infection")
          r <- data.table(ccs_epichron_label = colnames(get(l))[i],
                          n_positive_total = table(get(l)[sexo == s & agebands == age_label, ..i])[2],
                          n_positive_ctl = ifelse(length(get(l)[sexo == s & agebands == age_label, ..cols][, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 0, N]) == 0, 
                                                  yes = 0, 
                                                  no = get(l)[sexo == s & agebands == age_label, ..cols][, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 0, N]),
                          n_positive_hiv = ifelse(length(get(l)[sexo == s & agebands == age_label, ..cols][, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 1, N]) == 0, 
                                                  yes = 0, 
                                                  no = get(l)[sexo == s & agebands == age_label, ..cols][, .N, by = cols][get(cols[1]) == 1 & get(cols[2]) == 1, N]),
                          n_negative_total = table(get(l)[sexo == s & agebands == age_label, ..i])[1],
                          n_negative_ctl = ifelse(length(get(l)[sexo == s & agebands == age_label, ..cols][, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 0, N]) == 0, 
                                                  yes = 0, 
                                                  no = get(l)[sexo == s & agebands == age_label, ..cols][, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 0, N]),
                          n_negative_hiv = ifelse(length(get(l)[sexo == s & agebands == age_label, ..cols][, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 1, N]) == 0, 
                                                  yes = 0, 
                                                  no = get(l)[sexo == s & agebands == age_label, ..cols][, .N, by = cols][get(cols[1]) == 0 & get(cols[2]) == 1, N]),
                          # Odd ratio of the variable without intercept
                          OR = exp(coef(t))["HIV_infection1"],
                          # Confidence interval by Wald method
                          lb = ifelse(as.data.frame(exp(confint.default(t)))["2.5 %"]["HIV_infection1",] == Inf, NA, as.data.frame(exp(confint.default(t)))["2.5 %"]["HIV_infection1",]),
                          ub = ifelse(as.data.frame(exp(confint.default(t)))["97.5 %"]["HIV_infection1",] == Inf, NA, as.data.frame(exp(confint.default(t)))["97.5 %"]["HIV_infection1",]),
                          p_value = as.data.frame(summary(t)$coefficient)["Pr(>|z|)"]["HIV_infection1",]
          )
          r[is.na(ub) | is.na(lb), OR := NA]
          results <- rbind(results, r)
          
          },warning = function(w){message(paste0("no convergence at ", l," ", s, " ", colnames(get(l))[i]))}
        )

        rm(t)
      }
      # Remove empty row, round results and mark p values < 0.05:
      results <- results[!is.na(ccs_epichron_label)]
      
      cols <- grep("^n_positive_total|^OR|^lb|^ub|^p_value", colnames(results), value = T)
      results[is.na(n_positive_total), (cols) := NA]
      
      cols <- grep("^OR|^lb|^ub|^p_value", colnames(results), value = T)
      results[, (cols) := lapply(.SD, function(x) round(x, digits = 3)), .SDcols = cols]
      
      results[, p_05 := p_value < 0.05]
      
      results <- results[!(is.na(n_positive_total) | is.na(OR))]
      
      fwrite(results, paste0("output/",format(Sys.Date(), "%Y%m%d"),"/", l, "_", s, "_", gsub("<|>", "" ,age_label), "_", format(Sys.Date(), "%Y%m%d"), ".csv"))
    }
  }
}

# Save all results as one excel:
# Change to the result folder:
setwd(paste0("output/", format(Sys.Date(), "%Y%m%d")))
# Get file names
file.names = list.files(pattern="csv$", full.names = F)
# file.names <- file.names[-length(file.names)]
# Read them into a list
df.list = lapply(file.names, fread, encoding = "UTF-8")
names(df.list) <- gsub("diag_bool_|\\.csv$", "" ,file.names)

write.xlsx(df.list, paste0("Results_", format(Sys.Date(), "%Y%m%d"), ".xlsx"), overwrite = T)
setwd("../..")

