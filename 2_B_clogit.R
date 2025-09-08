# Author: Alejandro Santos Mejías
# Date last update: 2025-09-02
# Input: Boolean diagnosis datasets
#        Matched pair dataset
# Output: Conditional logistic regression output
# Motivational comment: Mistakes are unavoidable. 

rm(list = ls())
gc()

source("99_paths_and_packages.R")

library(survival)
library(dplyr)
library(purrr)

result_path <- paste0("output/", format(Sys.Date(),"%Y%m%d"))
if(!dir.exists(result_path)){dir.create(result_path)}

###########################################
##### Conditional Logistic regressions ####
###########################################

bdu <- fread("intermediate/bdu_full.csv", encoding = "UTF-8")
diag_bool_prev <- fread("intermediate/diagnosis_bool_prevalent_final.csv", encoding = "UTF-8")
cohort <- fread("intermediate/cohort_final.csv", encoding = "UTF-8")

# Assign match_id
t <- unique(cohort[patient_id %in% bdu$patient_id, .(patient_id)], by = "patient_id")
t[, match_id := 1:.N]
cohort <- merge(cohort, t, by = "patient_id")
fwrite(cohort,"intermediate/cohort_match_id.csv", encoding = "UTF-8")
cohort[, variable := NULL]
cohort <- melt(cohort, id.vars = c("vih_dt", "match_id"), value.name = "patient_id")
cohort <- unique(cohort, by = c("patient_id", "match_id"))

# Prepare variables for analysis in diagnosis dataset
diag_bool_prev <- merge(diag_bool_prev, bdu[, .(patient_id, edad, sexo, tsi, nacionalidad)], by = "patient_id")
age_cuts <- c(-Inf,44,65, Inf)
age_labels <- c("under 45", "45 to 65", "over 65")
diag_bool_prev[, agebands := cut(edad, age_cuts, labels = age_labels)]
diag_bool_prev[, edad := NULL]
setcolorder(diag_bool_prev, c("patient_id", "vih_dt", "HIV_infection", "sexo", "agebands", "tsi", "nacionalidad"))
diag_bool_prev[, nacionalidad := ifelse(nacionalidad == "ESPAÑA", "España", "Extranjero")]
diag_bool_prev[, nacionalidad := as.factor(nacionalidad)]
diag_bool_prev[, HIV_infection := as.factor(HIV_infection)]
diag_bool_prev[, tsi := factor(tsi, levels = c("< 18000", "entre 18000 y 100000", "> 100000", "Farmacia gratuita", "No asegurados", "Mutualistas"))]
setkey(diag_bool_prev, NULL)

# Duplicate controls in diagnosis dataset
# cohort[, .N, by = "patient_id"] # Check which one
diag_dup <- merge(diag_bool_prev, cohort[, .(patient_id, match_id)], by = "patient_id")
setcolorder(diag_dup, "match_id")

# Prepare formula for clogit models
# All comorbidities are alphabetically order after "nacionalidad" column, extract colnames of them
exposure_var <- colnames(diag_dup)[9:ncol(diag_dup)]
# Create all formulas for modelling
model_formula <- paste0("`", exposure_var,"`", " ~ HIV_infection + nacionalidad + tsi + strata(match_id)")

# Stratified analysis, so strata must be created in accordance with the previous created variables!!!
stratum <- setDT(expand.grid(sexo = c("HOMBRE", "MUJER"), edad = c("under 45", "45 to 65" , "over 65")))

for (i in 1:nrow(stratum)) {
  
  diag_subset <- copy(diag_dup[sexo == stratum[i, sexo] & agebands == stratum[i, edad]])
  
  model_results <- model_formula %>% map(.f = ~ clogit(formula = as.formula(.x), data = diag_subset)) %>%
    map(.f = ~ broom::tidy(.x, exponentiate = T, conf.int = T)) %>%
    bind_rows() %>%
    filter(term == "HIV_infection1") %>%
    mutate(across(where(is.numeric), round, digits = 3)) %>%
    mutate(ccs_epichron_label = exposure_var) %>%
    select(!c(term))
  
  setDT(model_results)
  setnames(model_results, c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high") ,c("OR", "std_error", "z_statistic", "p_value", "conf_lower", "conf_upper"))
  setcolorder(model_results, "ccs_epichron_label")
  
  # Create frequency table (this step can be transform to functional programming for sure)
  freq_results <- data.table(ccs_epichron_label = exposure_var, 
                             n_positive_total = diag_subset[, sapply(.SD, sum),  .SDcols = exposure_var],
                             n_positive_ctl = diag_subset[HIV_infection == 0, sapply(.SD, sum),  .SDcols = exposure_var],
                             n_positive_hiv = diag_subset[HIV_infection == 1, sapply(.SD, sum),  .SDcols = exposure_var])
  
  freq_results[, `:=`(
    n_negative_total = diag_subset[,.N] - n_positive_total,
    n_negative_ctl = diag_subset[HIV_infection == 0, .N] - n_positive_ctl,
    n_negative_hiv = diag_subset[HIV_infection == 1, .N] - n_positive_hiv
  )]
  
  # Generate final result table and prettify it
  results <- merge(freq_results, model_results, by = "ccs_epichron_label")
  results[conf_upper == Inf | is.na(conf_upper)| OR >= 1000, `:=`(OR = NA,
                                                      std_error = NA,
                                                      z_statistic = NA,
                                                      p_value = NA,
                                                      conf_lower = NA,
                                                      conf_upper = NA)]
  
  results[, p_05 := p_value < 0.05]
  results[, p_05 := as.character(p_05)][is.na(p_05), p_05 := "No convergence"]
    
  fwrite(results, paste0(result_path,"/clogit_full_", stratum[i, sexo], "_", gsub("<|>| ", "" , stratum[i, edad]), "_", format(Sys.Date(), "%Y%m%d"), ".csv"))
  
}

# Save all results as one excel:
# Change to the result folder:
setwd(result_path)
# Get file names
file.names = list.files(pattern="csv$", full.names = F)
# file.names <- file.names[-length(file.names)]
# Read them into a list
df.list = lapply(file.names, fread, encoding = "UTF-8")
names(df.list) <- gsub("^clogit_full_|_[0-9]*\\.csv$", "" ,file.names)

write.xlsx(df.list, paste0("Results_", format(Sys.Date(), "%Y%m%d"), ".xlsx"), overwrite = T)
setwd("../..")



