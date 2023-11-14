### plot_cloneGroRates.R, v1 | 14/11/19 ###
# Takes associated .RData objects and plots individual growth curves plots
#   for each strain; 
# Also plots combined max ODs, max rs, and lag phases for each strain.
# Modified from plot_mutantGroRates_v0.R

## Aesthetics to-do:
#

###

# load libraries & theme
library(ggplot2)
library(ggthemr)
ggthemr("pale")

# load data
load("data_cloneGroRates.RData")
# allData: all strains, including blanks
# ancData: clean plotting subset, only ancestral strains
# cloneData: as ancData, but only clone strains
# blankData: as ancdata, but only blank strains
# ancData_growths: list of summarized statistics for Ancestor max r[1], max OD[2], 
#   and lag-phase[3]; [[index]] accordingly to plot.

### GROWTH CURVES ###

# aesthetics
dAesBrk <- c("S0-R0", "S16-R0", "S0-R8", "S16-R8")
dAesLbl <- c("No drug", "STP", "RIF", "STP & RIF")
dAesClr <- c("S0-R0"="#18b518", "S16-R0"="#3e3bff","S0-R8"="#ff3e3b", "S16-R8"="#9f3d9d")
dAesGde <- guide_legend(title = "Environment")

# Plot 9 ancestors as a facetted plot
plot_ancestors <- ggplot(data=ancData, aes(x=time_hrs, y=blankedOD, col=drugs, group=interaction(drugs,rep))) +
  geom_line(size = 0.4) +
  scale_color_manual(values = dAesClr, breaks = dAesBrk, labels = dAesLbl, guide = dAesGde) +
  facet_wrap(~strain, ncol=2, dir="v") +
  ylim(0, 0.8) +
  scale_x_continuous(name = "Time (hours)", limits = c(0, 25)) +
  scale_y_continuous(name = "Optical Density", limits = c(0, 0.8)) +
  theme(axis.line = element_line(size = 0.5, colour = swatch()[1]),
        legend.position = c(0.78, 0.08),
        panel.border = element_rect(fill = NA, size = 1))
#plot_ancestors
plotName <- "ancestors"
ggsave(paste0(plotName, ".pdf"), plot = plot_ancestors, height=287, width=200, 
       units="mm", dpi = 300, limitsize = TRUE, path = "./plots")
ggsave(paste0(plotName, ".png"), plot = plot_ancestors, height=287, width=200, 
       units="mm", dpi = 600, limitsize = FALSE, path = "./plots")

# Plot all treatments, with 10 clones per plot
for(plate in unique(cloneData$plateID)){
  currentSubset <- cloneData[ cloneData$plateID == plate,]
  plot_clones <- ggplot(data=currentSubset, aes(x=time_hrs, y=blankedOD, col=drugs, group=interaction(drugs,rep))) +
    geom_line(size = 0.4) +
    scale_color_manual(values = dAesClr, breaks = dAesBrk, labels = dAesLbl, guide = dAesGde) +
    facet_wrap(~strain, ncol=2) +
    ylim(0, 0.8) +
    scale_x_continuous(name = "Time (hours)", limits = c(0, 25)) +
    scale_y_continuous(name = "Optical Density", limits = c(0, 0.8)) +
    theme(axis.line = element_line(size = 0.5, colour = swatch()[1]),
          panel.border = element_rect(fill = NA, size = 1))
  plotName <- plate
  ggsave(paste0(plotName, ".pdf"), plot = plot_clones, height=287, width=200, 
         units="mm", dpi = 300, limitsize = TRUE, path = "./plots")
  ggsave(paste0(plotName, ".png"), plot = plot_clones, height=287, width=200, 
         units="mm", dpi = 600, limitsize = FALSE, path = "./plots")
}

# # Plot all blanks in one plot, by columns
# plot_blanks <- ggplot(data=blankData, aes(x=time_hrs, y=OD, col=drugs, group=interaction(drugs,rep))) +
#   geom_line(size = 0.4) +
#   scale_color_manual(values = dAesClr, breaks = dAesBrk, labels = dAesLbl, guide = dAesGde) +
#   facet_wrap(~column, ncol=2) +
#   #ylim(0, 0.2) +
#   scale_x_continuous(name = "Time (hours)", limits = c(0, 25)) +
#   scale_y_continuous(name = "Optical Density", limits = c(0.085, 0.105)) +
#   theme(axis.line = element_line(size = 0.5, colour = swatch()[1]),
#         #legend.position = c(0.78, 0.08),
#         panel.border = element_rect(fill = NA, size = 1))
# plot_blanks
# plotName <- "blanks"
# ggsave(paste0(plotName, ".pdf"), plot = plot_blanks, height=287, width=200,
#        units="mm", dpi = 300, limitsize = TRUE, path = "./plots")
# ggsave(paste0(plotName, ".png"), plot = plot_blanks, height=287, width=200,
#        units="mm", dpi = 600, limitsize = FALSE, path = "./plots")

### MAX R, MAX OD, LAG-PHASE ###
# currently (13/11/19) only works for Ancestor plots
# ensure drugs and strains are refactored in the wrangler

# define ggplot aesthetics
dotBreaks <- c("IH_R24", "IH_R22", "IH_R27", "OH_RIF4",  "OH_RIF2", "OH_AB13RIF1", "HN_RIF3", "HN_RIF18", "HN_RIF7", "HN_RIF2", "OH_AB3STP1","IH_S07", "AB14", "HN_STP24", "HN_STP11", "OH_AB13STP3","AB3", "AB13")# don't use allStrains; this is to force S, R order
dotLabels <- c("IH_R24", "IH_R22", "IH_R27", "OH_RIF4",  "OH_RIF2", "OH_AB13RIF1", "HN_RIF3", "HN_RIF18", "HN_RIF7", "HN_RIF2", "OH_AB3STP1","IH_S07", "AB14", "HN_STP24", "HN_STP11", "OH_AB13STP3","AB3", "AB13")
dotColors <- c("IH_R24"="#660000", "IH_R22"="#FF0000", "IH_R27"="#FF9999", "OH_RIF4"="#994C00", "OH_RIF2"="#FF8000", 
               "OH_AB13RIF1"="#FFCC99", "HN_RIF3"="#999900", "HN_RIF18"="#FFFF00", "HN_RIF7"="#336600", "HN_RIF7"="#66CC00", 
               "OH_AB3STP1"="#33FF99", "IH_S07"="#99FFFF", "AB14"="#99FFCC", "HN_STP24"="#00CC00", "HN_STP11"="#99FF33", 
               "OH_AB13STP3"="#00CC66", "AB3"="#0080FF", "AB13"="#9999FF")
dotGuide <- guide_legend(title = "Strains")

# plot maximum growth rates
plot_maxr <- ggplot(data=ancData_growths[[1]], aes(x=strain, y=mean, col=strain)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0.15, size = 0.4, # "_bars" width = 0.5
                position = position_dodge(width = 0.75)) +
  #geom_vline(xintercept = c(0.5, 1.5, 2.5, 3.5, 4.5), color = "white") +
  scale_x_discrete(name = element_blank(), breaks = dotBreak_x, labels = dotLabel_x) +
  scale_y_continuous(name = expression("Maximum growth rate ["*min^-1*"]"), limits = c(0.0, 0.02), breaks = seq(0.000, 0.02, 0.005)) +
  #                     breaks = c(0.000, 0.0025, 0.005, 0.0075, 0.010, 0.0125)) +
  scale_color_manual(values = dotColors, breaks = dotBreaks, 
                     labels = dotLabels, guide = dotGuide) +
  facet_wrap(~drugs, ncol=4, dir="v") + 
  theme(axis.line = element_line(size = 0.5, colour = swatch()[1]),
        panel.border = element_rect(fill = NA, size = 1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
#plot_maxr
plotName <- "MaxRs"
ggsave(paste0(plotName, ".pdf"), plot = plot_maxr, dpi = 300, 
       height = 6, width = 8, limitsize = FALSE, path = "./plots")
ggsave(paste0(plotName, ".png"), plot = plot_maxr, dpi = 600, 
       height = 6, width = 8, limitsize = FALSE, path = "./plots")

# plot maxOD
plot_maxOD <- ggplot(data=ancData_growths[[2]], aes(x=strain, y=mean, col=strain)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.15, size = 0.4, 
                position = position_dodge(width = 0.75)) +
  #geom_vline(xintercept = c(0.5, 1.5, 2.5, 3.5, 4.5), color = "white") +
  scale_x_discrete(name = element_blank(), breaks = dotBreak_x, labels = dotLabel_x) +
  scale_y_continuous(name = expression("Maximum OD"[600]), limits = c(0.0, 0.8)) +
  scale_color_manual(values = dotColors, breaks = dotBreaks, 
                     labels = dotLabels, guide = dotGuide) +
  facet_wrap(~drugs, ncol=4, dir="v") + 
  theme(axis.line = element_line(size = 0.5, colour = swatch()[1]),
        panel.border = element_rect(fill = NA, size = 1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
#plot_maxOD
plotName <- "MaxODs"
ggsave(paste0(plotName, ".pdf"), plot = plot_maxOD, dpi = 300, 
       height = 6, width = 8, limitsize = FALSE, path = "./plots")
ggsave(paste0(plotName, ".png"), plot = plot_maxOD, dpi = 600, 
       height = 6, width = 8, limitsize = FALSE, path = "./plots")

# plot lag phase duration
plot_lagphase <- ggplot(data=ancData_growths[[3]], aes(x=strain, y=mean, col=strain)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.15, size = 0.4, 
                position = position_dodge(width = 0.75)) +
  #geom_vline(xintercept = c(0.5, 1.5, 2.5, 3.5, 4.5), color = "white") +
  scale_x_discrete(name = element_blank(), breaks = dotBreak_x, labels = dotLabel_x) +
  scale_y_continuous(name = "Lag Phase Duration [min]") +
  scale_color_manual(values = dotColors, breaks = dotBreaks, 
                     labels = dotLabels, guide = dotGuide) +
  facet_wrap(~drugs, ncol=4, dir="v") + 
  theme(axis.line = element_line(size = 0.5, colour = swatch()[1]),
        panel.border = element_rect(fill = NA, size = 1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
#plot_lagphase
plotName <- "LagPhases"
ggsave(paste0(plotName, ".pdf"), plot = plot_maxOD, dpi = 300, 
       height = 6, width = 8, limitsize = FALSE, path = "./plots")
ggsave(paste0(plotName, ".png"), plot = plot_maxOD, dpi = 600, 
       height = 6, width = 8, limitsize = FALSE, path = "./plots")


# end. #