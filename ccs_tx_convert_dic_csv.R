library(data.table)
library(stringr)

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
setnames(ccs, c("ccs_code", "ccs_label", "icd9_code"))
fwrite(ccs, "ccs_icd9_dic.csv")
