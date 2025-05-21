library(gtsummary)
library(data.table)
library(stringr)
library(ggplot2)
library(openxlsx)
library(psych) # Required for EFA
library(GPArotation) # Required for factor rotations

setwd("C:/Santos/VIH/Prueba_push/data7/")
options(scipen = 999)


##############
# Acknoledge #
##############
# For concepts: https://m-clark.github.io/posts/2020-04-10-psych-explained/#introduction
# For workflow: https://www.promptcloud.com/blog/exploratory-factor-analysis-in-r/
#               https://rpubs.com/pjmurphy/758265
# Original free source of the psych package author: https://www.personality-project.org/r/book/


######################################
##### Exploratory Factor Analysis ####
######################################

################################################################################

# Load the data
diag_bool_inc <- fread("intermediate/diagnosis_bool_incident_final.csv", encoding = "UTF-8")
diag_bool_prev <- fread("intermediate/diagnosis_bool_prevalent_final.csv", encoding = "UTF-8")
bdu <- fread("intermediate/bdu_full.csv", encoding = "UTF-8")
cohort <- fread("intermediate/cohort_final.csv", encoding = "UTF-8", colClasses = "character")

# Create parameters for filtering EFA:
prevalenceProportionVec <- 1:5
correlationKmo <- 0.5

# Create output folder:
if(!dir.exists(paste0("output/", format(Sys.Date(),"%Y%m%d")))){dir.create(paste0("output/", format(Sys.Date(),"%Y%m%d")))}
if(!dir.exists(paste0("output/", format(Sys.Date(),"%Y%m%d"), "/EFA_results"))){dir.create(paste0("output/", format(Sys.Date(),"%Y%m%d"), "/EFA_results"))}

for (prevalenceProportion in prevalenceProportionVec) {
  
# Filter out none-EFA-appropiate variables, variables with low number of
# positive cases or with really low correlations with the others.
# Remove variables with less than 1 positive case:
cols <- colnames(diag_bool_prev)[-c(1:3)]

for (i in 0:1) {
  t <- cols[diag_bool_prev[HIV_infection == i, colSums(.SD), .SDcols = cols]>2]
  assign(paste0("cols",ifelse(i == 0, "Ctl", "Hiv")), t)
  
}


# Test idoneity of data for EFA Kaiser-Meyer-Olkin factor adequacy is a measure
# of the proportion of variance among variables that might be common variance. A
# higher proportion indicates a higher KMO-value, which means the data is more
# suited for factor analysis. Rule of Thumbs:
# * 0.00 to 0.49 unacceptable
# * 0.50 to 0.59 miserable
# * 0.60 to 0.69 mediocre
# * 0.70 to 0.79 middling
# * 0.80 to 0.89 meritorious
# * 0.90 to 1.00 marvelous
for (i in 0:1) {
  
  kmo <- KMO(r = diag_bool_prev[HIV_infection == i, get(ifelse(i == 0, "colsCtl", "colsHiv")), with = FALSE])  
  if(i == 0){
    t <- colsCtl[kmo$MSAi > correlationKmo]
  }else{
    t <- colsHiv[kmo$MSAi > correlationKmo]
  }
  
  assign(paste0("cols",ifelse(i == 0, "Ctl", "Hiv")), t)
  
}

# Remove diseases with low incidence/prevalence. In theory for trustful result
# at least 10 obs per variable is desired, variables with less than 5 it is
# recommendable to remove. For prevalence, I will follow previous works,
# selecting prevalences > 1%, if one diseases index increase prevalence to up to
# 5%. For incidences, I will select diseases with at least 20 observations.
for (i in 0:1) {
  
  cols <- ifelse(i == 0, "colsCtl", "colsHiv")
  nObsDis <- diag_bool_prev[HIV_infection == i, .(disease = get(cols), N = colSums(.SD)/.N*100), .SDcols = get(cols)]
  t <- nObsDis[N >= prevalenceProportion, c(disease)] 
  assign(paste0("cols",ifelse(i == 0, "Ctl", "Hiv")), t)
}


# Let's start the analysis:
for (i in 0:1) {
  t <- ifelse(i == 0, "colsCtl", "colsHiv")
  cols <- get(t)

  # There are two possible workflows to elucidate MM patterns for a target disease
  # population:
  #  - 2 Factor analysis, one for the population with and without the target disease.
  #  - 1 Factor analysis for the whole cohort, checking in which factor the loading is greater for the target disease.
  
  # The first approach is followed in this analysis.
  
  # @param fm - factoring method, maximum likelihood chosen for consistency with previous work. minres the default one.
  # @param fa - fa = show eigen values for principal axis factor analysis.
  # @param cor - which correlation should be calculated, in our case tetrachoric (binary data).
  
  nFactors <-
    suppressMessages(suppressWarnings(fa.parallel(
      x = diag_bool_prev[HIV_infection == i, ..cols],
      fm = "ml",
      fa = "fa",
      cor = "tet"
    )))
  
  nMaxFactors <- nFactors$nfact
  
  if(nMaxFactors <= 2){
    warning("Optimal estimated factors less or equal 2. Check data!")
    nMaxFactors <- 3
  }
  
  # @param rotate - method to rotate the data for factoring. Take into account:
  # Orthogonal rotations if we assume our factors are NOT correlated. Such as none, varimax, quartimax, bentlerT, equamax, varimin, geominT and bifactor
  # Oblique rotations if we assume our factors are correlated. Such as oblimin, geominQ, Promax or target.rot
  for(l in 2:nMaxFactors){
   
    faResult <-
      suppressMessages(suppressWarnings(
        fa(
          diag_bool_prev[HIV_infection == i, ..cols],
          fm = "ml",
          nfactors = l,
          cor = "tet",
          rotate = "oblimin"
        )
      ))
    
    faLoadings <- as.data.frame(round(faResult$loadings[1:dim(faResult$loadings)[1], ], 3))
    faLoadings[faLoadings < 0.3] <- 0
    setorderv(faLoadings, colnames(faLoadings)[-1], rep(-1, length(colnames(faLoadings)[-1])), na.last = T)
    
    fwrite(faLoadings, paste0("output/", format(Sys.Date(),"%Y%m%d"), "/EFA_results/efa_loadings_", gsub("cols", "", t), "_", prevalenceProportion,"%_", l, "factors.csv"), encoding = "UTF-8", row.names = T)
    saveRDS(faResult, paste0("output/", format(Sys.Date(),"%Y%m%d"), "/EFA_results/efa_models_", gsub("cols", "", t), "_", prevalenceProportion,"%_", l, "factors.rds"))
     
  }
  
}
}
# TODO reliability statistic such as alpha and omega to check results

