library(data.table)
library(stringr)
library(ggplot2)

options(scipen = 999)

append_file <- function(directory = NA, pattern = NA, sep = "|", encoding = "UTF-8", label = T, colClasses = "character"){
  tables <- list.files(directory, pattern = pattern, full.names = T)
  result <- data.table()
  for (i in tables) {
    t <- fread(file = i, sep = sep, encoding = encoding, colClasses = colClasses)
    if(label){
      t[, ageband := str_extract(i, c("45to65|under45|over65"))]
      t[, sex := str_extract(i, c("HOMBRE|MUJER"))]
      }
    result <- rbindlist(list(result, t), fill = T)
  }
  return(result)
}

data <- append_file(directory = ".", pattern = "^clogit", sep = ",")

cols <- colnames(data)
cols <- cols[which(cols == "n_positive_total"):which(cols == "conf_upper")]
data[, (cols) := lapply(.SD, function(x) ifelse(x == "", NA, x)), .SDcols = cols]
data[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
data[, ageband := factor(ageband, levels = c("under45", "45to65", "over65"))]
data[, sex := factor(sex, levels = c("MUJER","HOMBRE"))]
data[, ccs_epichron_label := gsub("_", " ", ccs_epichron_label)]
# Define dictionary for labels
# Definimos el diccionario de traducciones
mapeo <- c(
  "Acute cerebrovascular disease" = "CVA",
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

# Aplicamos el cambio de forma masiva
data[, ccs_epichron_label := ifelse(ccs_epichron_label %in% names(mapeo), mapeo[ccs_epichron_label], ccs_epichron_label)]
setorder(data, ccs_epichron_label)

data[, ccs_epichron_label := factor(ccs_epichron_label, labels =  ,levels = unique(rev(ccs_epichron_label)))]
data[, p_05 := factor(p_05, levels = c("TRUE", "FALSE", "No convergence"))]

labels_full <- c(data[p_05 == "TRUE" & sex == "HOMBRE", unique(as.character(ccs_epichron_label))],
       data[p_05 == "TRUE" & sex == "MUJER", unique(as.character(ccs_epichron_label))])

labels_dup <- labels_full[duplicated(labels_full)]

labels_men <- data[p_05 == "TRUE" & sex == "HOMBRE", unique(as.character(ccs_epichron_label))]

labels_women <- data[p_05 == "TRUE" & sex == "MUJER", unique(as.character(ccs_epichron_label))]

labels_prev <- fread("significantPrev.csv", encoding = "UTF-8")
labels_prev <- labels_prev[, ccs_epichron_label := ifelse(ccs_epichron_label %in% names(mapeo), mapeo[ccs_epichron_label], ccs_epichron_label)][, ccs_epichron_label]

data2 <- copy(data[ccs_epichron_label %in% labels_dup[labels_dup %in% labels_prev]])

data2[ OR > 12, OR := 12]

data2[is.na(OR), c("OR", "conf_lower", "conf_upper") := 1.000001]

data2[, point := ifelse(data2$OR == 1.000001, 4, data2$ageband)]
data2[, point := factor(point, labels = c("under 45", "45 to 65", "over 65", "No convergence"))]
# Comment following line if ageband coloring:
data2[, point := factor(point, levels = c("No convergence","under 45", "45 to 65", "over 65"))]
# data2[, point := ifelse(data2$OR == 1.000001, "No convergence", "Convergence")]
# data2[, point := factor(point, levels = c("No convergence", "Convergence"))]

# t <- as.character(data2[ageband == "45to65"][order(OR), unique(ccs_epichron_label)])
# data2[, ccs_epichron_label := forcats::fct_relevel(data2$ccs_epichron_label, t)]

size = 1.5
dodge <- 0.5
# 
# #------------------------------------------------------------------------------#
# #### Facet by sex common OR in logistic and prevalences: ####
# #------------------------------------------------------------------------------#
# 
# p <- ggplot(data = data2, aes(x = OR, y = ccs_epichron_label, colour = ageband)) +
#   geom_vline(xintercept=1, lty=2) +
#   geom_vline(xintercept=c(0.1, 0.2, 0.5, 2, 5, 10), lty=3, colour = "gray70") +
#   # Finally working with pointrange being more esplicit
#   geom_pointrange(
#     aes(xmin = conf_lower, xmax = conf_upper, shape = point),
#     position = position_dodge(width = dodge),
#     linewidth = (size * 0.5)
#   ) +
#   labs(colour = "Conditions", shape = "Conditions") +
#   scale_y_discrete(
#     labels = function(x) {stringr::str_trunc(x, 25)},
#     expand = c(0, 1)
#   ) +
#   scale_shape_manual(values = c(
#     "No convergence" = 8,
#     "under 45" = 15,
#     "45 to 65" = 16,
#     "over 65" = 17
#   )) +
#   # Colour of the paper, switch to new journals
#   scale_colour_manual(values = c(
#     "under45" = "#ef3c40",
#     "45to65" = "#1cb7fa",
#     "over65" = "#308780"
#   )) +
#   # # Colour for sex:
#   # scale_colour_manual(values = c(
#   #   "HOMBRE" = "#ef3c40",
#   #   "MUJER" = "#1cb7fa"
#   # )) +
#   scale_x_continuous(
#     oob = scales::squish,
#     limits = c(0.1, 10),
#     breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
#     trans = "log10", # Transform to log10 scale, so all breaks are equally distanced
#     expand = c(0.05,0)
#   ) +
#   theme(
#     axis.text = element_text(color = "black"),
#     # text = element_text(size = 6),
#     axis.text.y = element_text(margin = margin(0, 0, 1, 0)),
#     strip.text.x = element_text(face = "bold", size = 10),
#     strip.text.y = element_text(face = "bold", size = 10),
#     # panel.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank(),
#     axis.title = element_blank(),
#     legend.position = "none",
#     # aspect.ratio = 2
#   ) +
#   # facet_grid(ageband~ sex)
#   facet_grid(~sex)
#   # facet_grid(~ageband)
# 
# p
# 
# # Max width = 19.05 (2250 px), width align with the text column 13.2, min width 6.68 cm (789 px)
# # Max height 22.23 cm, 2625 px
# ggsave(plot =  p , filename = "Prev_common_prueba2.tiff", width = 2250, height = 2625, unit = "px", dpi = 300)
# 
# 
# #------------------------------------------------------------------------------#
# #### Facet by sex common OR: ####
# #------------------------------------------------------------------------------#
# 
# data2 <- copy(data[ccs_epichron_label %in% labels_dup])
# 
# data2[ OR > 12, OR := 12]
# 
# data2[is.na(OR), c("OR", "conf_lower", "conf_upper") := 1.000001]
# 
# data2[, point := ifelse(data2$OR == 1.000001, 4, data2$ageband)]
# data2[, point := factor(point, labels = c("under 45", "45 to 65", "over 65", "No convergence"))]
# 
# p <- ggplot(data = data2, aes(x = OR, y = ccs_epichron_label, colour = ageband)) +
#   geom_vline(xintercept=1, lty=2) +
#   geom_vline(xintercept=c(0.1, 0.2, 0.5, 2, 5, 10), lty=3, colour = "gray70") +
#   # Finally working with pointrange being more esplicit
#   geom_pointrange(
#     aes(xmin = conf_lower, xmax = conf_upper, shape = point),
#     position = position_dodge(width = dodge),
#     linewidth = (size * 0.5)
#   ) +
#   labs(colour = "Conditions", shape = "Conditions") +
#   scale_y_discrete(
#     labels = function(x) {stringr::str_trunc(x, 25)},
#     expand = c(0, 1)
#   ) +
#   scale_shape_manual(values = c(
#     "No convergence" = 8,
#     "under 45" = 15,
#     "45 to 65" = 16,
#     "over 65" = 17
#   )) +
#   # Colour of the paper, switch to new journals
#   scale_colour_manual(values = c(
#     "under45" = "#ef3c40",
#     "45to65" = "#1cb7fa",
#     "over65" = "#308780"
#   )) +
#   scale_x_continuous(
#     oob = scales::squish,
#     limits = c(0.1, 10),
#     breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
#     trans = "log10", # Transform to log10 scale, so all breaks are equally distanced
#     expand = c(0.05,0)
#   ) +
#   theme(
#     axis.text = element_text(color = "black"),
#     # text = element_text(size = 6),
#     axis.text.y = element_text(margin = margin(0, 0, 1, 0)),
#     strip.text.x = element_text(face = "bold", size = 10),
#     strip.text.y = element_text(face = "bold", size = 10),
#     # panel.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank(),
#     axis.title = element_blank(),
#     legend.position = "none",
#     # aspect.ratio = 2
#   ) +
#   # facet_grid(ageband~ sex)
#   # facet_grid(~sex)
#   facet_grid(~ageband)
# 
# p
# 
# # Max width = 19.05 (2250 px), width align with the text column 13.2, min width 6.68 cm (789 px)
# # Max height 22.23 cm, 2625 px
# ggsave(plot =  p , filename = "all_prueba1.tiff", width = 2250, height = 2625, unit = "px", dpi = 300)

#------------------------------------------------------------------------------#
#### Facet by ageband common OR: ####
#------------------------------------------------------------------------------#

p <- ggplot(data = data2, aes(x = OR, y = ccs_epichron_label, colour = sex)) +
  geom_vline(xintercept=1, lty=2) +
  geom_vline(xintercept=c(0.1, 0.2, 0.5, 2, 5, 10), lty=3, colour = "gray70") +
  # Finally working with pointrange being more esplicit
  geom_pointrange(
    aes(xmin = conf_lower, xmax = conf_upper, shape = point),
    position = position_dodge(width = dodge),
    linewidth = (size * 0.5)
  ) +
  labs(colour = "Conditions", shape = "Conditions") +
  scale_y_discrete(
    labels = function(x) {stringr::str_trunc(x, 25)},
    expand = c(0, 1)
  ) +
  scale_shape_manual(values = c(
    "No convergence" = 8,
    "under 45" = 15,
    "45 to 65" = 16,
    "over 65" = 17
  )) +
  # # Colour of the paper, switch to new journals
  scale_colour_manual(values = c(
    "HOMBRE" = "#ef3c40",
    "MUJER" = "#1cb7fa"
  )) +
  scale_x_continuous(
    oob = scales::squish,
    limits = c(0.1, 10),
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
    trans = "log10", # Transform to log10 scale, so all breaks are equally distanced
    expand = c(0.05,0)
  ) +
  theme(
    axis.text = element_text(color = "black"),
    # text = element_text(size = 6),
    axis.text.y = element_text(margin = margin(0, 0, 1, 0)),
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 10),
    # panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    # aspect.ratio = 2
  ) +
  # facet_grid(ageband~ sex)
  # facet_grid(~sex)
  facet_grid( ~ ageband, labeller = as_labeller(c("under45" = "Under 45",
                                                  "45to65" = "Between 45 to 65",
                                                  "over65" = "Over 65")
  )
  )

p

# Max width = 19.05 (2250 px), width align with the text column 13.2, min width 6.68 cm (789 px)
# Max height 22.23 cm, 2625 px
ggsave(plot =  p , filename = "commonOR.png", width = 2250, height = 2625, unit = "px", dpi = 300)
ggsave(plot =  p , filename = "commonOR.tiff", width = 2250, height = 2625, unit = "px", dpi = 300)

# 
# #------------------------------------------------------------------------------#
# #### common OR by sex: ####
# #------------------------------------------------------------------------------#
# for (i in c("HOMBRE", "MUJER")) {
#   data2 <- copy(data[ccs_epichron_label %in% labels_dup & sex == i])
#   data2[OR > 12, OR := 12]
#   data2[is.na(OR), c("OR", "conf_lower", "conf_upper") := 1.000001]
#   data2[, point := ifelse(data2$OR == 1.000001, 4, data2$ageband)]
#   data2[, point := factor(point,
#                           labels = c("under 45", "45 to 65", "over 65", "No convergence"))]
#   
#   p <- ggplot(data = data2, aes(x = OR, y = ccs_epichron_label, colour = ageband)) +
#     geom_vline(xintercept = 1, lty = 2) +
#     geom_vline(
#       xintercept = c(0.1, 0.2, 0.5, 2, 5, 10),
#       lty = 3,
#       colour = "gray70"
#     ) +
#     # Finally working with pointrange being more esplicit
#     geom_pointrange(
#       aes(
#         xmin = conf_lower,
#         xmax = conf_upper,
#         shape = point
#       ),
#       position = position_dodge(width = dodge),
#       linewidth = (size * 0.5)
#     ) +
#     labs(colour = "Conditions", shape = "Conditions") +
#     scale_y_discrete(
#       labels = function(x) {
#         stringr::str_trunc(x, 25)
#       },
#       expand = c(0, 1)
#     ) +
#     scale_shape_manual(values = c(
#       "No convergence" = 8,
#       "under 45" = 15,
#       "45 to 65" = 16,
#       "over 65" = 17
#     )) +
#     # # Colour of the paper, switch to new journals
#     scale_colour_manual(values = c(
#       "under45" = "#ef3c40",
#       "45to65" = "#1cb7fa",
#       "over65" = "#308780"
#     )) +
#     scale_x_continuous(
#       oob = scales::squish,
#       limits = c(0.1, 10),
#       breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
#       trans = "log10",
#       # Transform to log10 scale, so all breaks are equally distanced
#       expand = c(0.05, 0)
#     ) +
#     theme(
#       axis.text = element_text(color = "black"),
#       # text = element_text(size = 6),
#       axis.text.y = element_text(margin = margin(0, 0, 1, 0)),
#       strip.text.x = element_text(face = "bold", size = 10),
#       strip.text.y = element_text(face = "bold", size = 10),
#       # panel.background = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.grid.major.x = element_blank(),
#       axis.title = element_blank(),
#       legend.position = "none",
#       # aspect.ratio = 2
#     ) +
#     # facet_grid(ageband~ sex)
#     # facet_grid(~sex)
#     facet_grid( ~ ageband, labeller = as_labeller(
#       c(
#         "under45" = "Under 45",
#         "45to65" = "Between 45 to 65",
#         "over65" = "Over 65"
#       )
#     ))
#   
#   p
#   
#   # Max width = 19.05 (2250 px), width align with the text column 13.2, min width 6.68 cm (789 px)
#   # Max height 22.23 cm, 2625 px
#   ggsave(
#     plot =  p ,
#     filename = paste0("all_prueba_", i, ".tiff"),
#     width = 2250,
#     height = 2625,
#     unit = "px",
#     dpi = 300
#   )
# }
# 

#------------------------------------------------------------------------------#
#### non-common OR HOMBRE: ####
#------------------------------------------------------------------------------#

data2 <- copy(data[ccs_epichron_label %in% labels_men[!(labels_men %in% labels_dup)] & sex == "HOMBRE"])
data2[, ageband := factor(ageband, levels = c("under45", "45to65", "over65"))]
data2[OR > 12, OR := 12]
data2[is.na(OR), c("OR", "conf_lower", "conf_upper") := 1.000001]
data2[, point := ifelse(data2$OR == 1.000001, 4, data2$ageband)]
data2[, point := factor(point,
                        labels = c("under 45", "45 to 65", "over 65", "No convergence"))]

p <- ggplot(data = data2, aes(x = OR, y = ccs_epichron_label, colour = ageband)) +
  geom_vline(xintercept=1, lty=2) +
  geom_vline(xintercept=c(0.1, 0.2, 0.5, 2, 5, 10), lty=3, colour = "gray70") +
  # Finally working with pointrange being more esplicit
  geom_pointrange(
    aes(xmin = conf_lower, xmax = conf_upper, shape = point),
    position = position_dodge(width = dodge),
    linewidth = (size * 0.5)
  ) +
  labs(colour = "Conditions", shape = "Conditions") +
  scale_y_discrete(
    labels = function(x) {stringr::str_trunc(x, 25)},
    expand = c(0, 1)
  ) +
  scale_shape_manual(values = c(
    "No convergence" = 8,
    "under 45" = 15,
    "45 to 65" = 16,
    "over 65" = 17
  )) +
  # # Colour of the paper, switch to new journals
  scale_colour_manual(values = c(
    "under45" = "#ef3c40",
    "45to65" = "#ef3c40",
    "over65" = "#ef3c40"
  )) +
  scale_x_continuous(
    oob = scales::squish,
    limits = c(0.1, 10),
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
    trans = "log10", # Transform to log10 scale, so all breaks are equally distanced
    expand = c(0.05,0)
  ) +
  theme(
    axis.text = element_text(color = "black"),
    # text = element_text(size = 6),
    axis.text.y = element_text(margin = margin(0, 0, 1, 0)),
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 10),
    # panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    # aspect.ratio = 2
  ) +
  # facet_grid(ageband~ sex)
  # facet_grid(~sex)
  facet_grid( ~ ageband, labeller = as_labeller(c("under45" = "Under 45",
                                                  "45to65" = "Between 45 to 65",
                                                  "over65" = "Over 65")
                                                )
              )

p

# Max width = 19.05 (2250 px), width align with the text column 13.2, min width 6.68 cm (789 px)
# Max height 22.23 cm, 2625 px
ggsave(plot =  p , filename = "hombreOR.tiff", width = 2250, height = 2625, unit = "px", dpi = 300)
ggsave(plot =  p , filename = "hombreOR.png", width = 2250, height = 2625, unit = "px", dpi = 300)

#------------------------------------------------------------------------------#
#### non-common OR MUJER: ####
#------------------------------------------------------------------------------#

data2 <- copy(data[ccs_epichron_label %in% labels_women[!(labels_women %in% labels_dup)] & sex == "MUJER"])
data2[, ageband := factor(ageband, levels = c("under45", "45to65", "over65"))]
data2[OR > 12, OR := 12]
data2[is.na(OR), c("OR", "conf_lower", "conf_upper") := 1.000001]
data2[, point := ifelse(data2$OR == 1.000001, 4, data2$ageband)]
data2[, point := factor(point,
                        labels = c("under 45", "45 to 65", "over 65", "No convergence"))]
data2 <- data2[ccs_epichron_label != "Cryptorchidism and other genitourinary congenital anomalies"]

p <- ggplot(data = data2, aes(x = OR, y = ccs_epichron_label, colour = ageband)) +
  geom_vline(xintercept=1, lty=2) +
  geom_vline(xintercept=c(0.1, 0.2, 0.5, 2, 5, 10), lty=3, colour = "gray70") +
  # Finally working with pointrange being more esplicit
  geom_pointrange(
    aes(xmin = conf_lower, xmax = conf_upper, shape = point),
    position = position_dodge(width = dodge),
    linewidth = (size * 0.5)
  ) +
  labs(colour = "Conditions", shape = "Conditions") +
  scale_y_discrete(
    labels = function(x) {stringr::str_trunc(x, 30)},
    expand = c(0, 1)
  ) +
  scale_shape_manual(values = c(
    "No convergence" = 8,
    "under 45" = 15,
    "45 to 65" = 16,
    "over 65" = 17
  )) +
  # # Colour of the paper, switch to new journals
  scale_colour_manual(values = c(
    "under45" = "#1cb7fa",
    "45to65" = "#1cb7fa",
    "over65" = "#1cb7fa"
  )) +
  scale_x_continuous(
    oob = scales::squish,
    limits = c(0.1, 10),
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
    trans = "log10", # Transform to log10 scale, so all breaks are equally distanced
    expand = c(0.05,0)
  ) +
  theme(
    axis.text = element_text(color = "black"),
    # text = element_text(size = 6),
    axis.text.y = element_text(margin = margin(0, 0, 1, 0)),
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 10),
    # panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    # aspect.ratio = 2
  ) +
  # facet_grid(ageband~ sex)
  # facet_grid(~sex)
  facet_grid( ~ ageband, labeller = as_labeller(c("under45" = "Under 45",
                                                  "45to65" = "Between 45 to 65",
                                                  "over65" = "Over 65")
  )
  )

p

# Max width = 19.05 (2250 px), width align with the text column 13.2, min width 6.68 cm (789 px)
# Max height 22.23 cm, 2625 px
ggsave(plot =  p , filename = "mujerORv2.tiff", width = 2250, height = 1200, unit = "px", dpi = 300)
ggsave(plot =  p , filename = "mujerORv2.png", width = 2250, height = 1200, unit = "px", dpi = 300)
