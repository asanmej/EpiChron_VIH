# Author: Alejandro Santos Mejías
# Date last update: 2025-04-25
# Input: Boolean diagnosis datasets
# Output: Top 20 incident and prevalent dieases
# Motivational comment: Life is a path full of suffering, you decide what to do with all of it! 


rm(list = ls())
gc()

source("99_paths_and_packages.R")

if(!dir.exists(paste0("output/", format(Sys.Date(),"%Y%m%d")))){dir.create(paste0("output/", format(Sys.Date(),"%Y%m%d")), recursive = TRUE)}

################################################################TRUE
##### Descriptive analysis #####
################################

################################################################################

bdu <- fread("intermediate/bdu_full.csv", encoding = "UTF-8")
diag_bool_prev <- fread("intermediate/diagnosis_bool_prevalent_final.csv", encoding = "UTF-8")
diag_bool_inc <- fread("intermediate/diagnosis_bool_incident_final.csv", encoding = "UTF-8")

# Let's prepare the sociodemographic database
bdu[, prevalente := bdu$patient_id %in% diag_bool_prev$patient_id]
bdu[, incidente := bdu$patient_id %in% diag_bool_inc$patient_id]
bdu[, nacionalidad := ifelse(nacionalidad == "ESPAÑA", "Spain", "Other")]

cols <- grep("dt$",colnames(bdu), value = T)
bdu[, (cols) := lapply(.SD, as.Date), .SDcols = cols]

cols <- grep("zbs_tipo|sexo",colnames(bdu), value = T)
bdu[, (cols) := lapply(.SD, as.factor), .SDcols = cols]

bdu[, tsi := factor(tsi, levels =  c("< 18000", "entre 18000 y 100000", "> 100000", "Farmacia gratuita", "Mutualistas", "No asegurados"), labels = c("< 18000 €/year ", "between 18000 to 100000 €/year", "> 100000 €/year", "Free medication", "Mutualist", "No health coverage"))]

bdu[, ageband := cut(edad, c(-Inf,44,65,Inf), labels = c("< 45 years", "45 - 65 years" , "> 65 years"))]

setnames(bdu, "indice_privación", "indice_privacion")

setcolorder(diag_bool_prev, c("patient_id", "vih_dt", "HIV_infection"))
setcolorder(diag_bool_inc, c("patient_id", "vih_dt", "HIV_infection"))

diag_bool_prev[, MM_prev := rowSums(diag_bool_prev[, -c(1:3)])]
diag_bool_inc[, MM_inc := rowSums(diag_bool_inc[, -c(1:3)])]

# BDU increases as there are multiple control patients repeated in disease dataset
bdu <- merge(bdu, diag_bool_prev[, .(patient_id, vih_dt, MM_prev)], by = "patient_id", all.x = T)
bdu <- merge(bdu, diag_bool_inc[, .(patient_id, vih_dt, MM_inc)], by = c("patient_id", "vih_dt"), all.x = T)
diag_bool_inc[, MM_inc := NULL]
diag_bool_prev[, MM_prev := NULL]

bdu <- bdu[!(prevalente == F & incidente == F)]
bdu[, muerte := ifelse(causa_baja == "FALLECIMIENTO", T, F)]


# Descriptive sociodemographic tables:
# To reduce code lines create a function, tweak the function to the needs.
theme_gtsummary_language("en", big.mark = "'")

descriptiveAnalysis <- function(diag, label) {
  diag %>%
    select(vih_bool,
           # nacimiento_dt,
           sexo,
           ageband,
           zbs_tipo,
           nacionalidad,
           tsi,
           MM_prev,
           muerte,
           indice_privacion) %>%
    mutate(
      vih_bool = ifelse(vih_bool == T, "HIV", "Non-HIV"),
      sexo = ifelse(sexo == "HOMBRE", "Men", "Women")
    ) %>%
    gtsummary::tbl_summary(
      by = vih_bool,
      label = list(
        ageband ~ "Ageband",
        # nacimiento_dt ~ "Birth date",
        sexo ~ "Sex",
        muerte ~ "Death",
        MM_prev ~ "Multimorbidity burden",
        zbs_tipo ~ "Living area",
        nacionalidad ~ "Nationality",
        tsi ~ "Copayment category",
        indice_privacion ~ "Deprivation Index"
      ),
      statistic = list(all_continuous() ~ "{mean} ({sd})"),
      digits = list(all_categorical() ~ c(0, 1)),
      missing = "no"
    ) %>%
    add_p(test = indice_privacion ~ "t.test") %>%
    bold_labels() %>%
    bold_p() %>%
    as_gt() %>% gt::gtsave(path = paste0("output/", format(Sys.Date(),"%Y%m%d")),
                           filename = paste0("Demographic_", label, ".html"))
}

# descriptiveAnalysis(diag = bdu[incidente == T], label = "inc")
descriptiveAnalysis(bdu[prevalente == T], "prev" )

# Multimorbidity by groups
bdu[, MM_prev := as.numeric(MM_prev)]
bdu[ageband == "< 45 years", MM_young := MM_prev]
bdu[ageband == "45 - 65 years", MM_mid := MM_prev]
bdu[ageband == "> 65 years", MM_old := MM_prev]

tmp <- bdu %>% select(vih_bool, sexo, MM_young, MM_mid, MM_old) %>% mutate(
  vih_bool = ifelse(vih_bool == T, "HIV", "Non-HIV"),
  sexo = ifelse(sexo == "HOMBRE", "Men", "Women")
) %>% tbl_strata(
  strata = sexo,
  .tbl_fun = ~ .x %>% tbl_summary(
    by = vih_bool,
    missing = "no",
    statistic = list(all_continuous() ~ "{mean} ({sd})"),
    label = list(
      MM_young = "Under 45 years",
      MM_mid = "Between 45 to 65 years",
      MM_old = "Over 65 years"
    )
  )  %>%
    add_p() %>%
    add_q(method = "bonferroni") %>%
    bold_p(q = T) %>%
    modify_column_hide(columns = "p.value") %>%
    modify_header(q.value ~ "**Adj. p-value**") %>%
    modify_footnote_header(
      footnote = "Wilcoxon rank sum test. Bonferrroni correction for multiple testing",
      columns = q.value,
      replace = T
    )
  ,
  .header = "**{strata}**\nN = {n}"
)

write.xlsx(bdu[, .(Media = mean(MM_prev,na.rm = T), SD = sd(MM_prev, na.rm = T)), by = c("vih_bool", "ageband", "sexo")], paste0("output/", format(Sys.Date(),"%Y%m%d"), "/multimorbidity_by_groups.xlsx"))
tmp %>% as_gt() %>% gt::gtsave(path = paste0("output/", format(Sys.Date(),"%Y%m%d")),
                               filename = paste0("multimorbidity_by_groups.html"))
tmp %>% as_hux_xlsx(paste0("output/", format(Sys.Date(),"%Y%m%d"), "/multimorbidity_by_groups.xlsx"))


# Top 20 prevalent diseases:
# cols <- colnames(diag_bool_prev)
# cols <- cols[-c(1:3)]
# cols <- cols[!unlist(diag_bool_prev[, lapply(.SD, function(x){sum(x) == 0}), .SDcols = cols, by = HIV_infection][, lapply(.SD, any), .SDcols = cols])]
# cols <- c("HIV_infection", cols)
# diag_bool_prev <- diag_bool_prev[, ..cols]
# 
# top20 <- data.table(
#   ccs = cols,
#   n_positive = t(diag_bool_prev[, lapply(.SD, sum), .SDcols = cols, by = "`HIV_infection`"][,-1]),
#   N_total_vih = table(diag_bool_prev$HIV_infection)[2],
#   N_total_ctl = table(diag_bool_prev$HIV_infection)[1],
#   p_value = t(diag_bool_prev[, lapply(.SD, function(x) round(chisq.test(x, diag_bool_prev$HIV_infection)$p.value, 4)), .SDcols = cols])
# )
# 
# setnames(top20, c("n_positive.V1", "n_positive.V2", "p_value.V1"), c("n_positive_ctl", "n_positive_vih", "p_value"))
# 
# top20[, p_value_adjusted := round(p.adjust(p_value, method = "bonferroni"), 3)]
# top20[, prop_vih := round(n_positive_vih/N_total_vih*100, 2)]
# top20[, prop_ctl := round(n_positive_ctl/N_total_ctl*100, 2)]
# top20[, diff := round(abs(prop_vih- prop_ctl), 3)]
# top20[, sign := ifelse((prop_vih - prop_ctl) >= 0, 1, -1)]
# setcolorder(top20, c("ccs", grep("vih$", colnames(top20), value = T), grep("ctl$", colnames(top20), value = T)))
# setorder(top20, -prop_vih)
# fwrite(top20, paste0("output/", format(Sys.Date(),"%Y%m%d"),"/Tabla_prevalencias_ranking.csv"))

# Prevalences stratified by sex
diag_bool_prev <- merge(diag_bool_prev, bdu[, .(patient_id, sexo)], by = "patient_id")
diag_bool_prev[, `:=`(patient_id = NULL, vih_dt = NULL)]

t <- diag_bool_prev[HIV_infection == 1 & sexo == "HOMBRE", lapply(.SD, function(x){sum(as.numeric(x), na.rm = T)/.N*100})]

setcolorder(diag_bool_prev, names(sort(unlist(t), decreasing = TRUE)))

t <- diag_bool_prev[sexo == "HOMBRE"] %>% mutate(
  Endometriosis = NA,
  Hypertension_complicating_pregnancy__childbirth_and_the_puerperium = NA
) %>% select(!sexo) %>%
  mutate(HIV_infection = ifelse(HIV_infection == 1, "HIV", "Non-HIV")) %>%
  tbl_summary(by = HIV_infection,
              missing = "no",
              missing_stat = "", 
              digits = list(all_categorical() ~ c(0,2))) %>%
  add_p() %>%
  add_q(method = "bonferroni") %>%
  bold_p(q = T) %>%
  modify_column_hide(columns = "p.value") %>%
  modify_header(q.value ~ "**Adj. p-value**") %>%
  modify_footnote_header(footnote = "Pearson’s Chi-squared test; Fisher’s exact test. Bonferrroni correction for multiple testing",
                         columns = q.value,
                         replace = T)

t2 <- diag_bool_prev[sexo == "MUJER"] %>% select(!sexo) %>%
  mutate(HIV_infection = ifelse(HIV_infection == 1, "HIV", "Non-HIV")) %>%
  tbl_summary(by = HIV_infection,
              missing = "no",
              missing_stat = "", 
              digits = list(all_categorical() ~ c(0,2))) %>%
  add_p() %>%
  add_q(method = "bonferroni") %>%
  bold_p(q = T) %>%
  modify_column_hide(columns = "p.value") %>%
  modify_header(q.value ~ "**Adj. p-value**") %>%
  modify_footnote_header(footnote = "Pearson’s Chi-squared test; Fisher’s exact test. Bonferrroni correction for multiple testing",
                         columns = q.value,
                         replace = T)


tmp <- tbl_merge(tbls = list(t, t2), tab_spanner = c("**Men**", "**Women**"))

tmp %>% as_gt() %>% gt::gtsave(path = paste0("output/", format(Sys.Date(), "%Y%m%d")),
                               filename = paste0("Prevalences_sex.html"))
tmp %>% as_hux_xlsx(paste0("output/", format(Sys.Date(),"%Y%m%d"), "/Prevalences_sex.xlsx"))


tmp_filtrado <- tmp %>%
  modify_table_body(~ .x %>%
                      dplyr::filter(q.value_1 < 0.05 | q.value_2 < 0.05))

tmp_filtrado %>% as_gt() %>% gt::gtsave(
  path = paste0("output/", format(Sys.Date(), "%Y%m%d")),
  filename = paste0("Prevalences_sex_filter.html")
)

tmp_filtrado %>% as_hux_table() %>%
  huxtable::quick_xlsx(
    file = paste0("output/",format(Sys.Date(), "%Y%m%d"),
    "/Prevalences_sex_filter.xlsx"
  ))
