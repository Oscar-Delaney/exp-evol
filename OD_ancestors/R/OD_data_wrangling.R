# load libraries
library(openxlsx)
library(growthrates)
library(tidyverse)

# source functions:
source("./R/functions.R")

#############################################################################
### Create tidy data frames using raw & spec files                        ###
#############################################################################

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
strainInfo <- read.csv("data/AB_resistant_mutants.csv") %>%
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
save(joined, file = "output/nice_OD_data.RData")

#############################################################################
### Run some simple statistical tests on the data                         ###
#############################################################################

# build a linear model to predict maxOD
maxODmodel <- lm(maxOD ~ Ancestor + CDS * STP + CDS * RIF, data = joined)
summary(maxODmodel)

# build a linear model to predict maxGrowthRate
maxGrowthRatemodel <- lm(maxGrowthRate ~ Ancestor + CDS * STP + CDS * RIF, data = joined)
summary(maxGrowthRatemodel)
