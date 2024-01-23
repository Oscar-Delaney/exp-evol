# Functions for wrangling, analysing and plotting optical density (OD) data.
# Code by Olivia Jessop, Oscar Delaney & Jan Engelstaedter

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
