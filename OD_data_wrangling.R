# load libraries
library(openxlsx)
library(growthrates)
library(tidyverse)

# takes raw_/spec_ files and transforms them into tidy data
transformODfile <- function(fileName, timePoints, timeInterval = 5) {
  # find starting position of data table within Excel file:
  column2 <- read.xlsx(paste0("./data/raw_", fileName), cols = 2, skipEmptyRows = FALSE, colNames = FALSE)[[1]]
  column2[is.na(column2)] <- "NA"
  startingRow <- match("Time", column2)
  # loading Excel file containing the raw read data as exported by the plate reader:
  rawdf <- read.xlsx(paste0("./data/raw_", fileName), 
                     rows = startingRow:(startingRow+timePoints+2), cols = 2:99, 
                     skipEmptyRows = FALSE)
  # loading Excel file containing specifications for all wells:
  specsdf <- read.xlsx(paste0("./data/specs_",fileName), cols = 1:97, skipEmptyRows = FALSE)
  
  # generating new data frame:
  times <- seq(0, by = timeInterval, length.out = timePoints)
  finaldf <- expand.grid(time = times, spec_idx = 1:96)

  # Populate the specs columns
  for (col in names(specsdf)) {
    finaldf[[col]] <- rep(specsdf[[col]], each = timePoints)
  }

  # Add OD data and other columns
  finaldf$OD <- unlist(rawdf[, -c(1:2)])
  return(finaldf)
}

# produces growth params
getGrowthParameters <- function(times, ODs, h = 5, tmax = Inf) {
  fit <- tryCatch(fit_easylinear(times[times < tmax], ODs[times < tmax], h = h), error=function(e) NULL)
  if (is.null(fit)) {
    return(c(y0 = NA, y0_lm = NA, mumax = NA, lag = NA, r2 = NA, deviance = NA, max = max(ODs)))
  } else {
    return(c(coef(fit), rsquared(fit), deviance = deviance(fit), max = max(ODs)))
  }
}

getBlankedODs <- function(allData) {
  allData %>%
    group_by(PlateName, Rep, STP, RIF) %>%
    mutate(median_blank_OD = median(OD[Strain == "BLANK"])) %>%
    mutate(blankedOD = pmax(0, OD - median_blank_OD)) %>%
    pull(blankedOD)
}

### Create tidy data frames using raw & spec files

allData <- data.frame() # create the allData frame before populating
fileNames <- gsub("^raw_", "", list.files(path = "./data", pattern = "^raw_.*\\.xlsx$"))
for (fileName in fileNames) { # populate allData
  trafoData <- transformODfile(fileName = fileName, timePoints = 24 * 60 / 5 + 1)
  allData <- rbind(allData, trafoData)
}

# allData$time_hrs <- allData$time / 60 # add column for time in hours
# allData$OD <- as.numeric(allData$OD) # make OD numeric
allData$blankedOD <- getBlankedODs(allData) # add column with blanked ODs

summarizedData <- allData %>%
  group_by(Strain, Rep, STP, RIF) %>%
  summarize(
    growthParams = list(getGrowthParameters(time, blankedOD, h = 20, tmax = 300)),
    .groups = 'drop'
  ) %>%
  rowwise() %>%
  mutate(
    maxGrowthRate = growthParams["mumax"],
    lagPhase = growthParams["lag"],
    rsquared = growthParams["r2"],
    maxOD = growthParams["max"]
  ) %>%
  select(-growthParams)

# Import strain info from AB_resistant_mutants.csv
strainInfo <- read.csv("AB_resistant_mutants.csv") %>%
  filter(CDS %in% c("rpoB", "rpsL")) %>%
  select(Mutant_Oscar_renamed, CDS, Ancestor, Mutation.Name)

joined <- left_join(summarizedData, strainInfo, by = c("Strain" = "Mutant_Oscar_renamed")) %>%
  # remove Strain rows with EMPTY or BLANK
  filter(!Strain %in% c("EMPTY", "BLANK")) %>%
  # if Ancestor is NA, set it equal to Strain
  mutate(Ancestor = ifelse(is.na(Ancestor), Strain, Ancestor),
    CDS = ifelse(is.na(CDS), "None", CDS)) %>%
  # if Mutation.Name is NA, set it equal to "WT"
  mutate(Mutation.Name = ifelse(is.na(Mutation.Name), "WT", Mutation.Name)) %>%
  select(c("Ancestor", "CDS", "Mutation.Name", "Rep", "STP", "RIF", "maxOD", "maxGrowthRate"))

# save as .RData:
save(joined, file = "summary.RData")

# build a linear model to predict maxOD
maxODmodel <- lm(maxOD ~ Ancestor + CDS * STP + CDS * RIF, data = joined)
summary(maxODmodel)

# build a linear model to predict maxGrowthRate
maxGrowthRatemodel <- lm(maxGrowthRate ~ Ancestor + CDS * STP + CDS * RIF, data = joined)
summary(maxGrowthRatemodel)
