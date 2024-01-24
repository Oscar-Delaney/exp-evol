# Load all the required libraries
library(openxlsx)
library(tidyverse)
library(stringr)
library(gganimate)
library(gifski)
library(ggnewscale)
library(patchwork)

# run functions:
source("./R/functions.R")

###############################################################################
### Wrangle the data:                                                       ###
###############################################################################

# Wrangle data for both com+ (=AB3) and com- (=AB13) strains:
dailyODDataAB3 <- wrangle_dailies("./data/specs/spec_290720-070820_for_2023.xlsx",
                                  "./data/ODs/dailyreads AB3/")
dailyODDataAB13 <- wrangle_dailies("./data/specs/spec_290720-070820_for_2023.xlsx",
                                   "./data/ODs/dailyreads AB13/")

# merge and save the data:
dailyODData <- rbind(dailyODDataAB3 |> 
                       mutate(Strain = "AB3", Rec = "com+"),
                     dailyODDataAB13 |> 
                       mutate(Strain = "AB13", Rec = "com-")) |>
  select(Strain, Rec, everything())

write_csv(dailyODData, "./output/OD_data_wrangled.csv")

###############################################################################
### Produce plots:                                                          ###
###############################################################################

# plots over all transfers:
plotDailyODs(dailyODData,
             "./plots/DailyODs_All.pdf", height = 12, width = 7)

# plot only endpoints for individual treatments:
plotDailyODs(dailyODData |> 
               filter(Strain == "AB3", Treatment == "CMB", transfer == "T10"), 
             "./plots/EndpointODs_AB3_CMB.pdf", width = 3, height = 3)
plotDailyODs(dailyODData |> 
               filter(Strain == "AB3", Treatment == "MIX", transfer == "T10"), 
             "./plots/EndpointODs_AB3_MIX.pdf", width = 3, height = 3)
plotDailyODs(dailyODData |> 
               filter(Strain == "AB13", Treatment == "CMB", transfer == "T10"), 
             "./plots/EndpointODs_AB13_CMB.pdf", width = 3, height = 3)
plotDailyODs(dailyODData |> 
               filter(Strain == "AB13", Treatment == "MIX", transfer == "T10"), 
             "./plots/EndpointODs_AB13_MIX.pdf", width = 3, height = 3)
plotEndpointODs(dailyODData |> filter(Treatment != "NOD"),
                "./plots/EndpointODs.pdf", width = 8, height = 4)

# Animated plots for presentations:
animateDailyODs(dailyODDataAB3, "./plots/DailyODs_AB3.gif")
animateDailyODs(dailyODDataAB13, "./plots/DailyODs_AB13.gif")

# Combination & mixing treatment only:
animateDailyODs(dailyODDataAB3 |> filter(Treatment == "CMB"), "./plots/DailyODs_AB3_CMB.gif",
                width = 3, height = 3, units = "in", res = 200, pointsize = 6)
animateDailyODs(dailyODDataAB13 |> filter(Treatment == "CMB"), "./plots/DailyODs_AB13_CMB.gif",
                width = 3, height = 3, units = "in", res = 200, pointsize = 6)
animateDailyODs(dailyODDataAB3 |> filter(Treatment == "MIX"), "./plots/DailyODs_AB3_MIX.gif",
                width = 3, height = 3, units = "in", res = 200, pointsize = 6)
animateDailyODs(dailyODDataAB13 |> filter(Treatment == "MIX"), "./plots/DailyODs_AB13_MIX.gif",
                width = 3, height = 3, units = "in", res = 200, pointsize = 6)
