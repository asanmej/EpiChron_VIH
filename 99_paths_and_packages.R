library(gtsummary)
library(data.table)
library(lubridate)
library(stringr)
library(ggplot2)
library(openxlsx)
library(psych) # Required for EFA
library(GPArotation) # Required for factor rotations
library(plotly)
library(htmlwidgets)

setwd("C:/Santos/vih_Datos8")
options(scipen = 999)

source("mapeo_ciap_cie9/mapping_ciap_bifap_cie9_v3.0.R")