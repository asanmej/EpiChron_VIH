library(data.table)
library(stringr)
library(tidyr)

setwd("C:/Santos/VIH/datos7/ccs_r/")

ccs <- readLines("AppendixASingleDX.txt")
ccs <- str_wrap(ccs)
ccs <- str_replace_all(ccs, "\n", " ")

ccs[ccs == ""] <- "\t"
ccs[grepl("[a-zA-Z]$|\\)$", ccs)] <- paste(ccs[grepl("[a-zA-Z]$|\\)$", ccs)],"%%")

ccs <- str_flatten(ccs, collapse = "__")
ccs <- str_split(ccs, pattern = "\t__")[[1]]
ccs <- str_replace_all(ccs, "__", " ")
ccs <- str_replace_all(ccs, " $", "")
ccs <- str_split(ccs, " %% ")
ccs <- as.data.table(t(as.data.table(ccs)))
ccs <- ccs[-1]
setnames(ccs, c("ccs_label", "ICD9_codes"))
ccs <- data.table(str_split_fixed(ccs$ccs_label, "(?<=[0-9]) ", n = 2), ICD9_codes = ccs$ICD9_codes)
setnames(ccs, c("ccs_code", "ccs_label", "cie9_no_dot"))
ccs <- separate_longer_delim(ccs, c(cie9_no_dot), delim = " ")
# ccs[, icd9_code := paste0("^",str_replace_all(ccs$icd9_code, " ", "|^"))]

# Add multilevel ccs mapping:
mltLevel <- fread("ccs_multi_dx_tool_2015.csv", encoding = "UTF-8", quote = "'", colClasses = "character")
setnames(mltLevel, c("cie9_no_dot", paste0(c("ccs_code_lvl", "ccs_label_lvl"), rep(1:4, each = 2))))

cols <- colnames(mltLevel)
mltLevel[, (cols) := lapply(.SD, function(x) gsub("\\\"", "", x)), .SDcols = cols]
mltLevel[, (cols) := lapply(.SD, function(x) ifelse(x == " ", NA, x)), .SDcols = cols]
mltLevel[, cie9_no_dot := str_wrap(cie9_no_dot)]

t <- merge(ccs, mltLevel, by = "cie9_no_dot", all.x = T)
setDT(t)
t[, ccs_code := as.integer(ccs_code)]

# Add clinical updates by Luis and Jonas:
t[, ccs_epichron_code := ccs_code]
t[cie9_no_dot == "2729", ccs_epichron_code := 53]
t[cie9_no_dot == "71650", ccs_epichron_code := 203]
t[cie9_no_dot == "73730", ccs_epichron_code := 209]
t[cie9_no_dot == "27800", ccs_epichron_code := 300]
t[cie9_no_dot == "3319", ccs_epichron_code := 653]

t[cie9_no_dot %in% c("27549", "2759"), ccs_epichron_code := 301]
t[cie9_no_dot %in% c("30751","30759","3071"), ccs_epichron_code := 302]
t[cie9_no_dot %in% c("30272","30276","30275","30270"), ccs_epichron_code := 303]
t[cie9_no_dot %in% c("3068","30747","30746"), ccs_epichron_code := 304]
t[cie9_no_dot %in% c("3007", "30081"), ccs_epichron_code := 305]

t[ccs_code == 101, ccs_epichron_code := 100]

# Drop the following codes:
t <- t[!(ccs_code %in% c(124, 181, 663) | cie9_no_dot %in% c("3070", "2510", "2512", "6963", "6965"))]

# Redistribute some ccs_categories. In case of grouping, first code use as main. For segregation, a free code from >300 was chosen:
t[ccs_code %in% c(49, 50), ccs_epichron_code := 49] # Group together Diabetes
t[ccs_code %in% c(98, 99), ccs_epichron_code := 98] # Group together chronic hypertensions
t[ccs_code %in% 11:47, ccs_epichron_code := 11] # Group together all neoplasms
t[ccs_code_lvl3 == "3.11.2", ccs_epichron_code := 300] # Separate Obesity from the rest of other nutritional and metabolic disorders
setkey(t, NULL)

# Assign and create new labels for ccs_epichron:
temp <- unique(t[, .(ccs_code, ccs_label)])
setnames(temp, c("ccs_epichron_code", "ccs_epichron_label"))
t <- merge(t, temp, by = "ccs_epichron_code", all.x = T)
t[ccs_epichron_code == 49, ccs_epichron_label := "Diabetes"]
t[ccs_epichron_code == 98, ccs_epichron_label := "Hypertension"]
t[ccs_epichron_code == 11, ccs_epichron_label := "Neoplasms"]
t[ccs_epichron_code == 300, ccs_epichron_label := "Obesity"]

# Update some labels with clearer concepts:
t[ccs_epichron_code ==51, ccs_epichron_label := "Hyperparathyroidism and other endocrine disorders"]
t[ccs_epichron_code ==54, ccs_epichron_label := "Gout"]
t[ccs_epichron_code ==59, ccs_epichron_label := "Thalassemia and other anemia"]
t[ccs_epichron_code ==81, ccs_epichron_label := "Restless leg syndrome and other nervous system conditions"]
t[ccs_epichron_code ==89, ccs_epichron_label := "Vision defects"]
t[ccs_epichron_code ==90, ccs_epichron_label := "Escleritis"]
t[ccs_epichron_code ==91, ccs_epichron_label := "Myodesopsia and other eye disorders"]
t[ccs_epichron_code ==93, ccs_epichron_label := "Menieres disease"]
t[ccs_epichron_code ==94, ccs_epichron_label := "Hearing loss and other ear and sense organ disorders"]
t[ccs_epichron_code ==95, ccs_epichron_label := "Neuropathy and other nervous system disorders"]
t[ccs_epichron_code ==100, ccs_epichron_label := "Ischemic heart disease"]
t[ccs_epichron_code ==134, ccs_epichron_label := "Allergic rhinitis"]
t[ccs_epichron_code ==151, ccs_epichron_label := "Hepatic steatosis and other liver diseases"]
t[ccs_epichron_code ==155, ccs_epichron_label := "Irritable bowel and other GI disorders"]
t[ccs_epichron_code ==163, ccs_epichron_label := "Urinary incontinence"]
t[ccs_epichron_code ==175, ccs_epichron_label := "Dyspareunia"]
t[ccs_epichron_code ==198, ccs_epichron_label := "Psoriasis and other inflammatory"]
t[ccs_epichron_code ==211, ccs_epichron_label := "Polymyalgia rheumatica"]
t[ccs_epichron_code ==214, ccs_epichron_label := "Tongue tie an other digestive congenital anomalies"]
t[ccs_epichron_code ==215, ccs_epichron_label := "Cryptorchidism and other genitourinary congenital anomalies"]
t[ccs_epichron_code ==253, ccs_epichron_label := "Atopic dermatitis"]
t[ccs_epichron_code ==259, ccs_epichron_label := "Obstructive Sleep apnea"]
t[ccs_epichron_code ==301, ccs_epichron_label := "Other microcrystalline arthritis"]
t[ccs_epichron_code ==302, ccs_epichron_label := "Eating disorders"]
t[ccs_epichron_code ==303, ccs_epichron_label := "Sexual disorders"]
t[ccs_epichron_code ==304, ccs_epichron_label := "Sleeping disorders"]
t[ccs_epichron_code ==305, ccs_epichron_label := "Somatization and hypochondria disorders"]
t[ccs_epichron_code ==661, ccs_epichron_label := "Adictions"]

t[ccs_epichron_label == "Mood disorders", ccs_epichron_label := "Depression and Mood disorders"]

fwrite(t, "ccs_icd9_dic.csv")

