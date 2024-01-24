# load libraries
# install grow96 if necessary via:
# devtools::install_github("JanEngelstaedter/grow96")
library(grow96)
library(growthrates)
library(tidyverse)

# source functions:
source("./R/functions.R")

#############################################################################
### Load and process OD data using raw & spec files                       ###
#############################################################################

dat <- processODData(specPath = "./data", dataPath = "./data") |>
  blankODs()
# produce a quality control report:
qcODData(dat, path = "./output/qc")

# add meta data:
dat <- dat |>
  left_join(read_csv("./data/mutant_key.csv") |> select(Strain, com, CDS, Mutation_name), 
            by = join_by(Strain)) |>
  select(Plate:Strain, com, CDS, Mutation_name, everything()) |>
  mutate(Mutation_name = ifelse(Strain %in% c("AB3", "AB13"), 
                                "wildtype", 
                                Mutation_name)) |>
  mutate(com = ifelse(Strain %in% c("AB3", "AB13"), 
                      c(AB3 = "com+", AB13 = "com-")[Strain], 
                      com))
  

# calculate statistics such as maximum OD and growth rate:
datAnalysed <- analyseODData(dat, tmax = 12 * 60)

# correct maximum growth rate to NA when maxOD is too low:
datAnalysed <- lapply(datAnalysed, FUN = function(x) {x$mumax <- ifelse(x$maxOD < 0.08, NA, x$mumax)
                                                      return(x) })

#############################################################################
### Run some simple statistical tests on the data                         ###
#############################################################################

# build a linear model to predict maxOD
maxODmodel <- lm(maxOD ~ com * CDS * Drug, data = datAnalysed$pars)
summary(maxODmodel)

# build a linear model to predict maxGrowthRate
maxGrowthRatemodel <- lm(mumax ~ com * CDS * Drug, data = datAnalysed$pars)
summary(maxGrowthRatemodel)


#############################################################################
### Produce plots                                                         ###
#############################################################################

# define ggplot aesthetics:
mutations <- c("wildtype", "rpsL_A128C", "rpsL_A262G", "rpsL_G275C",
               "rpoB_T443G", "rpoB_T1555C", "rpoB_C1556T", "rpoB_A1568C", "rpoB_C1712T")
colors <- c(hsv(0.33, 0.7, 0.8), hsv(seq(0.6, 0.75, length.out = 3), 0.7, 1),
            hsv(seq(0.15, 0, length.out = 5), 0.7, 1))

# transform columns into factors for plotting:
datAnalysed$pars <- datAnalysed$pars %>%
  mutate(Mutation_name = factor(Mutation_name, levels = mutations)) %>%
  mutate(com = factor(com, levels = c("com+", "com-"))) %>%
  mutate(Drug = factor(Drug, levels = c("No drug", "STP", "RIF", "STP + RIF")))

datAnalysed$means <- datAnalysed$means %>%
  mutate(Mutation_name = factor(Mutation_name, levels = mutations)) %>%
  mutate(com = factor(com, levels = c("com+", "com-"))) %>%
  mutate(Drug = factor(Drug, levels = c("No drug", "STP", "RIF", "STP + RIF")))

# plot maximum OD:
p_maxOD <- ggplot() +
  geom_point(data = datAnalysed$pars,
             mapping = aes(x = Mutation_name, y=maxOD, col = Mutation_name, shape = com),
             position = position_dodge(width = 0.66), size = 0.4) +
  geom_point(data = datAnalysed$means,
             mapping = aes(x = Mutation_name, y=maxOD, col = Mutation_name, shape = com),
             position = position_dodge(width = 0.66)) +
  facet_wrap(~Drug, ncol=4, dir="v") +
  scale_x_discrete(name = element_blank(), labels = element_blank()) +
  scale_y_continuous(name = expression("Maximum optical density")) +
  scale_color_manual(values = colors, labels = mutations, name = "Mutation") +
  scale_shape_manual(values = c("com+" = 19, "com-" = 1), name = "Ancestor") +
  theme_bw()

ggsave("./plots/maxODs.pdf", plot = p_maxOD, height = 5, width = 8)


# plot maximum growth rates:
p_maxGrowth <- ggplot() +
  geom_point(data = datAnalysed$pars,
             mapping = aes(x = Mutation_name, y=mumax, col = Mutation_name, shape = com),
             position = position_dodge(width = 0.66), size = 0.4) +
  geom_point(data = datAnalysed$means,
             mapping = aes(x = Mutation_name, y=mumax, col = Mutation_name, shape = com),
             position = position_dodge(width = 0.66)) +
  facet_wrap(~Drug, ncol=4, dir="v") +
  scale_x_discrete(name = element_blank(), labels = element_blank()) +
  scale_y_continuous(name = expression("Maximum growth rate ["*h^-1*"]")) +
  scale_color_manual(values = colors, labels = mutations, name = "Mutation") +
  scale_shape_manual(values = c("com+" = 19, "com-" = 1), name = "Ancestor") +
  theme_bw()

ggsave("./plots/maxGrowthRates.pdf", plot = p_maxGrowth, height = 5, width = 8)

# plotting the raw growth curves:
p_growth <- ggplot(data = dat |> filter(WellType == "DATA") |>
         mutate(Mutation_name = factor(Mutation_name, levels = mutations)) %>%
         mutate(com = factor(com, levels = c("com+", "com-"))) %>%
         mutate(Drug = factor(Drug, levels = c("No drug", "STP", "RIF", "STP + RIF")))) +
  geom_line(aes(x = Time_h, y = blankedOD, col = com, group = interaction(Replicate, com))) +
  facet_grid(vars(Mutation_name), vars(Drug)) +
  scale_color_manual(values = c("com+" = "#1b9e77", "com-" = "#7570b3"), name = "Ancestor") +
  xlab("Time [h]") +
  ylab("OD600 (blanked)") +
  theme_bw() +
  theme(legend.position="top")

ggsave("./plots/growthCurves.pdf", plot = p_growth, height = 10, width = 8)
