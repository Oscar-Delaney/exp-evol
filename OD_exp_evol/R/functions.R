# Functions to process and plot OD data from the evolution experiment
# Code by Olivia Jessop, Michael Thompson & Jan Engelstaedter

#' Subtract the blanking weights from the OD data. Don't want blanking? just set to 0
#' @param odVec a vector containing all the original OD data
#' @param treatment a string of the treatment
#' @return a new vector with updated values after blanking
blanking <- function(odVec, treatment) {
  # Store all the blanking values
  # For now user to update them here
  STR = c(rep(0.085966667, 3), rep(0.08665, 3))
  RIF = c(rep(0.085966667, 3), rep(0.0899, 3))
  CMB = c(rep(0.085966667, 3), rep(0.09, 3))
  MIX = c(rep(0.085966667, 3), 0.0899, 0.08665, 0.09)
  NOD = c(rep(0.085966667, 3), rep(0.085966667, 3))
  
  # Create a data frame for string lookup
  df <- data.frame(STR = STR, RIF = RIF, CMB = CMB, MIX = MIX, NOD = NOD)
  
  # Retrieve the correct blanking values and return the updated vector
  blankingWeights <- df[treatment][[1]]
  return (odVec - blankingWeights)
}

#' Creates a dataframe containing daily data based on a spec file
#' @param specFile a string containing the location of the spec file
#' @param folderPath a string containing the location of the data
#' @return A dataframe containing all the wrangled data from the excel files
wrangle_dailies <- function(specFile, folderPath, blanking = FALSE) {
  
  #' Read in the data from the specfile
  datesGuide <- read.xlsx(specFile, sheet = 1, skipEmptyRows = FALSE)
  plateGuide <- read.xlsx(specFile, sheet = 2, skipEmptyRows = FALSE)
  sheetGuide <- read.xlsx(specFile, sheet = 3, skipEmptyRows = FALSE)
  
  # generate total plates from number of sheets in sheet guide
  plateCount <- length(sheetGuide$Sheet)
  
  # create main dataframe; loop outputs get appended here
  finaldf <- data.frame()
  
  # Add a progress bar so we know it's running
  print("Wrangling")
  pb <- txtProgressBar(min = 1, max = length(datesGuide$Date), initial = 1)
  
  # LOOP 1 - INDIVIDUAL DAILY FILES #
  for(i in 1:length(datesGuide$Date)){
    # Progress bar!
    setTxtProgressBar(pb,i)
    
    # specify current file name vars, then parse current filename
    currentDate <- datesGuide$Date[i] # date
    currentTrf <- datesGuide$TransferID[i] # T-number
    rawFileName = paste0(folderPath, "/", "dailyread_", currentTrf, "_", currentDate, ".xlsx")
    
    # create current data frame; this is appended to finaldf
    currentdf <- data.frame(date = rep(currentDate, (24*plateCount)), 
                            transfer = rep(currentTrf, (24*plateCount)), 
                            plateID = NA, wellID = NA, wellCol = NA, wellRow = NA,
                            Treatment = NA, Replicate = NA, Environment = NA, OD = NA)
    
    # LOOP 2 -  INDIVIDUAL PLATES WITHIN EACH FILE #
    for(plate in 1:plateCount){
      # offsets the entries for each loop for adding to the dataframe
      offsetVec <- ((plate-1) * 24 + 1):(plate*24)
      
      # Populate the data frame
      currentdf$plateID[offsetVec] <- rep(sheetGuide$PlateID[plate], 24)
      currentdf$wellID[offsetVec] <- plateGuide$Well[1:24]
      currentdf$wellCol[offsetVec] <- plateGuide$Column[1:24]
      currentdf$wellRow[offsetVec] <- plateGuide$Row[1:24]
      currentdf$Treatment[offsetVec] <- rep(sheetGuide$PlateType[plate],24)
      # old code: WRONG!!!
      #currentdf$Replicate[offsetVec] <- rep(1:4 + 4*(plate+1) %% 2, each = 6)
      # corrected version:
      # infer shift from last letter of plate name:
      shift <- as.integer(stringr::str_sub(sheetGuide$PlateID[plate], -1))
      currentdf$Replicate[offsetVec] <- rep(1:4 + (shift - 1) * 4, each = 6)
      
      odData <- as.vector(t(as.matrix(slice(read.xlsx(rawFileName, sheet = plate, rows = 3:7, cols = 3:8, colNames = TRUE, skipEmptyRows = FALSE), c(1,2,3,4)))))
      # Blanking
      if (blanking) {
        blankedOD <- blanking(odData, sheetGuide$PlateType[plate])
        currentdf$OD[offsetVec] <- blankedOD
      } else {
        currentdf$OD[offsetVec] <- odData
      }
      currentdf$Environment[((plate-1) * 24 + 1):(plate*24)] <- plateGuide$Environment[ plateGuide$PlateType == sheetGuide$PlateType[plate] ]
    }
    close(pb)
    # append currentdf to finaldf
    finaldf <- rbind(finaldf, currentdf)
    
  }
  fullTypes <- c("STR" = "STP only", "RIF" = "RIF only", "CMB" = "Combination", "MIX" = "Mixed", "NOD" = "No drug")
  finaldf$Treatment_f <- fullTypes[finaldf$Treatment]
  
  # factorize <Treatment> to force facet label order
  finaldf$Treatment_f <- factor(finaldf$Treatment_f, levels = as.vector(fullTypes))
  
  return(finaldf)
}

#' Plot all the data
#' @param allData a dataframe with all the wrangled data
plotDailyODs <- function(allData, fileName = NULL, ...) {
  if (allData |> pull(Strain) |> unique() |> length() == 1L) {
    # only one of either com+ or com- to plot ==> facet by treatment
    p <- ggplot(data=allData, aes(wellCol, Replicate, fill = OD)) +
      geom_tile(colour="black",size=0.20) +
      coord_equal() +
      scale_x_continuous(name = "Subpopulation", breaks = seq(1, 6, 1), position = "top", expand = expansion(0, 0)) +
      scale_y_reverse(name = "Population", breaks = seq(1, 8, 1), expand = expansion(0, 0)) + 
      scale_fill_continuous(high = "dodgerblue4", low = "white", limits = c(0.0, 2.0), breaks = seq(0.0, 2.0, 0.5)) +
      theme(legend.position = "bottom", panel.spacing = unit(0.1, "lines")) +
      facet_grid(Treatment_f~transfer, switch = "x")
  } else {
    # both com+ and com- to plot ==> facet by com and have different plots for each treatment
    treats <- allData |> pull(Treatment_f) |> unique()
    p <- NULL
    for(i in 1:length(treats)) {
      p <- p / ggplot(data=allData |> filter(Treatment_f == treats[i]), 
                      aes(wellCol, Replicate, fill = OD)) +
        geom_tile(colour="black",size=0.20) +
        coord_equal() +
        scale_x_continuous(name = NULL, breaks = seq(1, 6, 1), position = "top", expand = expansion(0, 0)) +
        scale_y_reverse(name = "Population", breaks = seq(1, 8, 1), expand = expansion(0, 0)) + 
        scale_fill_continuous(high = "dodgerblue4", low = "white", limits = c(0.0, 2.0), breaks = seq(0.0, 2.0, 0.5)) +
        theme(legend.position = "bottom", panel.spacing = unit(0.1, "lines")) +
        facet_grid(Rec~transfer, switch = "x") +
        ggtitle(treats[i])
    }
    p <- p +
      plot_layout(guides = 'collect', axis_titles = "collect") &
      theme(legend.position = "bottom")
  }
  
  plot(p)
  if (!is.null(fileName)) {
    ggsave(fileName, plot = p, dpi = 300, 
           ...,
           limitsize = FALSE)
  }
}

#' Plot endpoint data
#' @param allData a dataframe with all the wrangled data
plotEndpointODs <- function(allData, 
                            fileName = NULL, ...) {
  # merge and filter data:
  endpointData <- allData |>
    filter(transfer == "T10")
  
  # add dummy for subpopulation legend:
  legendData <- endpointData |>
    filter(Replicate == 1, Rec == "com-")
  legendData$OD <- 0
  legendData$Rec <- ""
  legendData$Replicate <- 0
  legendData$Environment <- factor(legendData$Environment,
                                   levels = c("X", "S", "R", "B"))
  
  allData <- rbind(endpointData, legendData)
  allData$Rec <- factor(allData$Rec, levels = c("", "com+", "com-"))
  
  # colours for the four environments:
  envColours <- c("#FEFAE8", "#D52D20", "#0095FF", "#68339A")
  names(envColours) <- c("X", "R", "S", "B")
  
  facetplot <- ggplot() +
    scale_x_discrete(name = "Subpopulation", position = "bottom", labels = NULL) +
    scale_y_reverse(name = "Replicate population", breaks = seq(1, 8, 1), expand = expansion(0, 0)) + 
    geom_tile(data = legendData,
              mapping = aes(as.factor(wellCol), Replicate, fill = Environment), 
              colour="black") +
    scale_fill_manual(values = envColours, labels = c("No drug", "STP", "RIF", "Both")) +
    new_scale_fill() +
    scale_fill_continuous(high = "dodgerblue4", low = "white", limits = c(0.0, 2.0), breaks = seq(0.0, 2.0, 0.5)) +
    geom_tile(data = endpointData,
              mapping = aes(as.factor(wellCol), Replicate, fill = OD), 
              colour="black", size=0.20) +
    theme(legend.position = "right", panel.spacing = unit(0.1, "lines")) +
    facet_grid(Rec ~ Treatment_f, scales = "free_y", space = "free_y") +
    theme_minimal()
  plot(facetplot)
  if (!is.null(fileName)) {
    ggsave(fileName, plot = facetplot, dpi = 300, 
           ...,
           limitsize = FALSE)
  }
}


#' Animate all the data into a GIF
#' @param allData a dataframe with all the wrangled data
animateDailyODs <- function(allData, fileName = NULL, ...) {
  allData$transfer <- as.numeric(gsub("T", "", allData$transfer))
  animPlot <- ggplot(data=allData, aes(wellCol, Replicate, fill = OD)) +
    geom_tile(colour="black", linewidth=0.2) +
    coord_equal() +
    scale_x_continuous(name = "Well", breaks = seq(1, 6, 1), position = "top", expand = expansion(0, 0)) +
    scale_y_reverse(name = "Population", breaks = seq(1, 8, 1), expand = expansion(0, 0)) + 
    scale_fill_continuous(high = "dodgerblue4", low = "white", limits = c(0.0, 2.0), breaks = seq(0.0, 2.0, 0.5)) +
    theme(legend.position = "bottom", panel.spacing = unit(0.1, "lines")) +
    facet_grid(~Treatment_f, switch = "x") +
    labs(title = 'Transfer: {closest_state}') +
    transition_states(transfer, state_length = 5, wrap = FALSE) +
    ease_aes('linear')
  anim <- animate(animPlot, renderer = gifski_renderer(loop = FALSE), ...) #, duration = 30, end_pause = 80, ...)
  anim_save(fileName, anim)
  return(anim)
}

