# Corrected script for forest plots of multiple logistic models
# Fixed empty syntax argument in factor conversion and improved str_extract implementation

library(data.table)
library(stringr)
library(gt)
library(dplyr)

options(scipen = 999)

# result_path <- paste0("output/", format(Sys.Date(),"%Y%m%d"))
result_path <- "output/20260714_art1_clogit_results"
if(!dir.exists(result_path)){dir.create(result_path)}

# Read and append all clogit result tables from the directory
append_file <- function(directory = ".", pattern = "^clogit", sep = ",", encoding = "UTF-8", label = TRUE, colClasses = "character") {
  tables <- list.files(directory, pattern = pattern, full.names = TRUE)
  result <- data.table()
  
  for (i in tables) {
    t <- fread(file = i, sep = sep, encoding = encoding, colClasses = colClasses)
    if (label) {
      # Correctly extract a single matched pattern from the file path
      t[, ageband := str_extract(i, "45to65|under45|over65")]
      t[, sex := str_extract(i, "HOMBRE|MUJER")]
    }
    result <- rbindlist(list(result, t), fill = TRUE)
  }
  return(result)
}

data <- append_file(directory = result_path, pattern = "^clogit", sep = ",")

# Convert numeric columns to appropriate types
cols <- colnames(data)
cols <- cols[which(cols == "n_positive_total"):which(cols == "conf_upper")]
data[, (cols) := lapply(.SD, function(x) ifelse(x == "" | is.na(x), NA, x)), .SDcols = cols]
data[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]

# Set factor variables
data[, ageband := factor(ageband, levels = c("under45", "45to65", "over65"))]
data[, sex := factor(sex, levels = c("MUJER", "HOMBRE"))]
data[, ccs_epichron_label := gsub("_", " ", ccs_epichron_label)]

# Define dictionary mapping for long labels
mapeo <- c(
  "Acute cerebrovascular disease" = "CVA",
  "Adictions" = "Addictions",
  "Attention-deficit  conduct  and disruptive behavior disorders" = "Attention-deficit",
  "Delirium  dementia  and amnestic and other cognitive disorders" = "Delirium dementia",
  "Diseases of white blood cells" = "Diseases of WBC",
  "Disorders of lipid metabolism" = "Dyslipemia",
  "Infective arthritis and osteomyelitis  except that caused by tuberculosis or sexually transmitted disease " = "Infective arthritis and others",
  "Inflammatory conditions of male genital organs" = "Inflammatory genital conditions",
  "Congestive heart failure  nonhypertensive" = "CHF",
  "Coagulation and hemorrhagic disorders" = "Coagul. and hmrrg diseases",
  "Chronic obstructive pulmonary disease and bronchiectasis" = "COPD and bronchiectasis", 
  "Depression and Mood disorders" = "Depression and others",
  "Headache  including migraine" = "Headache",
  "Hepatic steatosis and other liver diseases" = "Hepatic steatosis",
  "HIV-related disease" = "Infective HIV-related disease",
  "Hyperparathyroidism and other endocrine disorders" = "Hyperparathyroidism",
  "Neuropathy and other nervous system disorders" = "Neuropathy and others",
  "Other male genital disorders" = "Others genital disorders",
  "Other non-traumatic joint disorders" = "Non-traumatic join",
  "Peri-  endo-  and myocarditis  cardiomyopathy  except that caused by tuberculosis or sexually transmitted disease " = "Peri- and endo-myocarditis",
  "Peripheral and visceral atherosclerosis" = "Atherosclerosis",
  "Psoriasis and other inflammatory" = "Psoriasis",
  "Restless leg syndrome and other nervous system conditions" = "Restless leg and others",
  "Retinal detachments  defects  vascular occlusion  and retinopathy" = "Retinopathy",
  "Thalassemia and other anemia" = "Thalassemia",
  "Transient cerebral ischemia" = "TIA"
)

# Apply label mapping mass update
data[, ccs_epichron_label := ifelse(ccs_epichron_label %in% names(mapeo), mapeo[ccs_epichron_label], ccs_epichron_label)]
setorder(data, ccs_epichron_label)

# Fixed syntax: removed empty 'labels =' argument in factor declaration
data[, ccs_epichron_label := factor(ccs_epichron_label, levels = unique(rev(ccs_epichron_label)))]
data[, p_05_adj := factor(p_05_adj, levels = c("TRUE", "FALSE", "No convergence"))]

data <- data[p_05_adj == "TRUE"]

# Adjust formatting of OR and CI
data[, OR_CI := sprintf("%.2f (%.2f, %.2f)", OR, conf_lower, conf_upper)]
# Set proper p value labels
data[, p_val_formatted := "<0.001"]
# Additional adjustments:
data[, ageband := factor(ageband, levels = c("under45", "45to65", "over65"), labels = c("Under 45 years", "Between 45 to 65 years", "Over 65 years"))]
data[, n_total := n_positive_total + n_negative_total]
data[, n_hiv_total := n_positive_hiv + n_negative_hiv]
data[, ccs_epichron_label := as.character(ccs_epichron_label)]
setorder(data, -OR)

# Set up variable names
dt_formatted <- data[, .(
  sex,
  ageband,
  OR,
  Condition = ccs_epichron_label,
  `n cases` = n_positive_total,
  `Total N` = n_total,
  `n PLWH cases` = n_positive_hiv,
  `Total PLWH N` = n_hiv_total,
  `OR (95% CI)` = OR_CI,
  `Bonferroni p-value` = p_val_formatted
)]


build_stratified_table <- function(data_input, target_sex, title_text) {
  data_input %>%
    filter(sex == target_sex) %>%
    arrange(ageband, desc(OR)) %>%
    select(-sex, -OR) %>%
    gt(groupname_col = "ageband") %>%
    # tab_header(title = title_text) %>%
    cols_label(
      Condition = "Condition",
      `n cases` = "n cases",
      `Total N` = "Total N",
      `n PLWH cases` = "n PLWH cases",
      `Total PLWH N` = "Total PLWH N",
      `OR (95% CI)` = "OR (95% CI)",
      `Bonferroni p-value` = "Bonferroni p-value"
    ) %>%
    cols_align(
      align = "center",
      columns = c(`n cases`, `Total N`, `n PLWH cases`, `Total PLWH N`, `OR (95% CI)`, `Bonferroni p-value`)
    ) %>%
    cols_align(
      align = "left",
      columns = Condition
    ) %>%
    tab_options(
      row_group.background.color = "#EFEFEF",
      row_group.font.weight = "bold",
      table.font.size = px(12)
    )
}


table_men_gt <- build_stratified_table(dt_formatted, "HOMBRE", "Significant Conditions by Age Stratum (Males)")
table_men_gt

gtsave(table_men_gt, paste0(result_path, "/logit_result_table_men.html"))

table_women_gt <- build_stratified_table(dt_formatted, "MUJER", "Significant Conditions by Age Stratum (Females)")
table_women_gt
gtsave(table_women_gt, paste0(result_path, "/logit_result_table_women.html"))

# Forest plot included in the tables:

plot_data <- data[p_05_adj == "TRUE"]

generate_publication_forestplot <- function(data_input, target_sex, output_filename) {
  
  df <- data_input[sex == target_sex]
  setorder(df, ageband, -OR)
  
  dt_list <- list()
  age_levels <- levels(df$ageband)
  
  for (ag in age_levels) {
    sub_df <- df[ageband == ag]
    if (nrow(sub_df) > 0) {
      header_row <- data.table(
        Condition = as.character(ag),
        cases_text = "",
        plwh_text = "",
        or_ci_text = "",
        OR = NA_real_,
        conf_lower = NA_real_,
        conf_upper = NA_real_,
        is_header = TRUE
      )
      
      data_rows <- data.table(
        Condition = paste0("  ", sub_df$ccs_epichron_label),
        cases_text = sprintf("%d / %d", sub_df$n_positive_total, sub_df$n_total),
        plwh_text = sprintf("%d / %d", sub_df$n_positive_hiv, sub_df$n_hiv_total),
        or_ci_text = sprintf("%.2f (%.2f, %.2f)", sub_df$OR, sub_df$conf_lower, sub_df$conf_upper),
        OR = sub_df$OR,
        conf_lower = sub_df$conf_lower,
        conf_upper = sub_df$conf_upper,
        is_header = FALSE
      )
      
      dt_list[[ag]] <- rbind(header_row, data_rows)
    }
  }
  
  final_df <- rbindlist(dt_list)
  final_df$` ` <- paste(rep(" ", 20), collapse = " ")
  
  setnames(final_df, 
           old = c("Condition", "cases_text", "plwh_text", "or_ci_text"),
           new = c("Condition", "Total Cases\n(n / N)", "PLWH Cases\n(n / N)", "OR (95% CI)"))
  
  table_data <- final_df[, .(`Condition`, `Total Cases\n(n / N)`, `PLWH Cases\n(n / N)`, ` `, `OR (95% CI)`)]
  
  tm <- forest_theme(
    base_size = 8,
    ci_pch = 16,
    ci_col = "black",
    ci_fill = "black",
    ci_alpha = 1,
    ci_line_col = "black",
    refline_col = "black",
    refline_lty = "dashed",
    vertline_col = "black",
    vertline_lty = "dotted",
    row_gpar = gpar(fill = "white"),
    
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 9),
      bg_params = list(fill = "white")
    ),
    
    page_margin = unit(c(0, 0, 0, 0), "pt")
  )
  
  p <- forest(
    table_data,
    est = final_df$OR,
    lower = final_df$conf_lower,
    upper = final_df$conf_upper,
    sizes = 0.4,
    is_summary = final_df$is_header,
    ci_column = 4,
    ref_line = 1,
    vert_line = 2,
    x_trans = "log10",
    ticks_at = c(0.1, 0.5, 1, 2, 5, 10),
    xlim = c(0.1, max(final_df$conf_upper, na.rm = TRUE) * 1.2),
    theme = tm
  )
  
  # Center text and header of all columns but first one
  p <- edit_plot(p, col = 2:5, which = "text", hjust = unit(0.5, "npc"), x = unit(0.5, "npc"))
  p <- edit_plot(p, col = 2:5, part = "header", which = "text", hjust = unit(0.5, "npc"), x = unit(0.5, "npc"))
  
  header_indices <- which(final_df$is_header)
  for (idx in header_indices) {
    p <- edit_plot(p, row = idx, col = 1, part = "body", which = "text",
                   gp = gpar(fontface = "bold", fontsize = 9))
    
    for (col_idx in 1:5) {
      p <- edit_plot(p, row = idx, col = col_idx, part = "body", 
                     which = "background", gp = gpar(fill = "#EFEFEF"))
    }
    
    p <- add_border(p, row = idx, col = 1:5, where = "top", gp = gpar(col = "black", lwd = 1))
    p <- add_border(p, row = idx, col = 1:5, where = "bottom", gp = gpar(col = "black", lwd = 1))
  }
  
  dims <- get_wh(p, unit = "in")
  
  crop_factor_w <- 0.98
  crop_factor_h <- 0.93
  
  width_in <- dims[1] * crop_factor_w
  height_in <- dims[2] * crop_factor_h
  
  png(output_filename, width = width_in, height = height_in, units = "in", res = 300)
  grid.newpage()
  pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
  grid.draw(p)
  popViewport()
  dev.off()
  
  return(p)
}

fp_men <- generate_publication_forestplot(plot_data, "HOMBRE", paste0(result_path, "/forest_plot_table_men.png"))
fp_women <- generate_publication_forestplot(plot_data, "MUJER", paste0(result_path, "/forest_plot_table_women.png"))

plot(fp_men)
plot(fp_women)