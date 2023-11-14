# load libraries
library(openxlsx)
library(growthrates)
library(tidyverse)

### FUNCTIONS ###

# takes raw_/spec_ files and transforms them into tidy data
transformODfile <- function(fileName, timePoints, timeInterval = 5, 
                            saveAsRData = FALSE, saveAsXLSX = FALSE) {
  # find starting position of data table within Excel file:
  column2 <- read.xlsx(paste0("raw_", fileName), cols = 2, skipEmptyRows = FALSE, colNames = FALSE)[[1]]
  column2[is.na(column2)] <- "NA"
  startingRow <- match("Time", column2)
  # loading Excel file containing the raw read data as exported by the plate reader:
  rawdf <- read.xlsx(paste0("raw_", fileName), 
                     rows = startingRow:(startingRow+timePoints+2), cols = 2:99, 
                     skipEmptyRows = FALSE)
  # loading Excel file containing specifications for all wells:
  specsdf <- read.xlsx(paste0("specs_",fileName), cols = 1:97, skipEmptyRows = FALSE)
  
  # generating new data frame:
  times <- seq(0, (timePoints - 1) * timeInterval, by = timeInterval)
  finaldf <- data.frame(time = rep(times, 96), date = NA, rep = NA, plateID = NA, reader = NA, 
                        row = NA, column = NA, well = NA, strain = NA, STP = NA, RIF = NA, OD = NA)
  
  # filling the other columns of the data frame:
  for(i in 1:96) {
    finaldf$date[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$Date[i], timePoints)
    finaldf$rep[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$Rep[i], timePoints)
    finaldf$plateID[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$PlateName[i], timePoints)
    finaldf$reader[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$Reader[i], timePoints)
    finaldf$row[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$Row[i], timePoints)
    finaldf$column[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$Column[i], timePoints)
    finaldf$well[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$Well[i], timePoints)
    finaldf$strain[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$Strain[i], timePoints)
    finaldf$STP[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$STP[i], timePoints)
    finaldf$RIF[((i-1) * timePoints + 1):(i * timePoints)] <- rep(specsdf$RIF[i], timePoints)
    finaldf$OD[((i-1) * timePoints + 1):(i * timePoints)] <- rawdf[, i + 2]
  }
  if (saveAsRData == TRUE) save(finaldf, file = paste0("./transformed outputs/trafo_", fileName, ".RData"))
  if (saveAsXLSX == TRUE) write.xlsx(finaldf, file = paste0("./transformed outputs/trafo_", fileName, ".xlsx"))
  return(finaldf)
} 

# produces growth params
getGrowthParameters <- function(times, ODs, h = 5, tmax = Inf) {
  fit <- tryCatch(fit_easylinear(times[times<tmax], ODs[times<tmax], h = h), error=function(e) NA)
  if (is.na(fit)) return(c(y0 = NA, y0_lm = NA, mumax = NA, lag = NA, r2 = NA, deviance = NA, max = max(ODs)))
  else return(c(coef(fit), rsquared(fit), deviance = deviance(fit), max = max(ODs)))
} # tryCatch will return NAs instead of errors, ESPECIALLY if "growthrates" isn't loaded

# produce vector of blanked ODs, specific to plate, time point and drug environment:
getBlankedODs <- function(allData) {
  blankedODs <- rep(NA, nrow(allData))
  for (plateID in unique(allData$plateID)) {
    indices <- (allData$plateID == plateID)
    for(rep in unique(allData$rep)) {
      indices <- (indices & allData$rep == rep)
      for(drugs in unique(allData$drugs)) {
        indices <- (indices & allData$drugs == drugs)
        for(time in unique(allData$time)) {
          indices <- (indices & allData$time == time)
          meanBlank <- mean(allData$OD[indices & allData$strain == "BLANK"])
          blankedODs[indices] <- allData$OD[indices] - meanBlank
        }
      }
    }
  }
  return(blankedODs)
}


# takes summarized growths and summarizes again for stats, and refactors for plotting
# editing this is tedious, need to generalize
summarize_for_plotting <- function(data){
  
  # summarize a bunch of statistics for plotting and error bars
  summarizedGrowthRates <- summarize(group_by(data, strain, drugs), 
                                     mean = mean(maxGrowthRate),
                                     se = sd(maxGrowthRate)/sqrt(length(maxGrowthRate)))
  summarizedODs <- summarize(group_by(data, strain, drugs), 
                             mean = mean(maxOD),
                             se = sd(maxOD)/sqrt(length(maxOD)))
  summarizedLagPhases <- summarize(group_by(data, strain, drugs), 
                                   mean = mean(lagPhase),
                                   se = sd(lagPhase)/sqrt(length(lagPhase)))
  
  # refactor to force facet/point orders
  summarizedGrowthRates$strain <- factor(summarizedGrowthRates$strain,
                                         levels = c("IH_R24", "IH_R22", "IH_R27", "OH_RIF4",  "OH_RIF2", "OH_AB13RIF1", "HN_RIF3", "HN_RIF18", "HN_RIF7", "HN_RIF2", "OH_AB3STP1","IH_S07", "AB14", "HN_STP24", "HN_STP11", "OH_AB13STP3","AB3", "AB13")) 
  summarizedODs$strain <- factor(summarizedODs$strain,
                                 levels = c("IH_R24", "IH_R22", "IH_R27", "OH_RIF4",  "OH_RIF2", "OH_AB13RIF1", "HN_RIF3", "HN_RIF18", "HN_RIF7", "HN_RIF2", "OH_AB3STP1","IH_S07", "AB14", "HN_STP24", "HN_STP11", "OH_AB13STP3","AB3", "AB13"))
  summarizedLagPhases$strain <- factor(summarizedLagPhases$strain,
                                       levels = c("IH_R24", "IH_R22", "IH_R27", "OH_RIF4",  "OH_RIF2", "OH_AB13RIF1", "HN_RIF3", "HN_RIF18", "HN_RIF7", "HN_RIF2", "OH_AB3STP1","IH_S07", "AB14", "HN_STP24", "HN_STP11", "OH_AB13STP3","AB3", "AB13"))
  
  return(list(summarizedGrowthRates, summarizedODs, summarizedLagPhases))
}

### FUNCTIONS END ###
