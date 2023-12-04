# 04/12/2023
# megamodel, code for the multifunctionality manuscript
# N.A. Pichon

library(lme4)

# load All_BigData
# load plot_description_scaled, which includes SLA and MPD values
# source the multidiv function
# source the stepba function

# setwd("~/")
All_BigData <- read.table("All_BigData.txt", header = T)

source("stepba function.R")
source("multidiv.R")

#

##### scale and data trasformation ####


All_BigData[c("Species_richness", "Nitrogen", "Fungicide", 
              "CWM_SLA", "MPD_SLA_abundance", "MPD_SLA_presence")] = scale(All_BigData[c("Species_richness", "Nitrogen", "Fungicide", 
                                                                                         "CWM_SLA", "MPD_SLA_abundance", "MPD_SLA_presence")])

All_BigData$Aboveground_biomass = sqrt(All_BigData$Aboveground_biomass)
All_BigData$Herbivory = sqrt(All_BigData$Herbivory)
# Pathogens no transformation
All_BigData$Soil_respiration = log(All_BigData$Soil_respiration)
All_BigData$Plant_N_uptake = log(All_BigData$Plant_N_uptake)
All_BigData$Plant_P_uptake = log(All_BigData$Plant_P_uptake)
All_BigData$Belowground_biomass = log(All_BigData$Belowground_biomass)
All_BigData$BGlucosidase = sqrt(All_BigData$BGlucosidase)
All_BigData$Phosphatase = log(All_BigData$Phosphatase)
All_BigData$Carbon_storage = sqrt(All_BigData$Carbon_storage+1)


##### MegaData ####

# Calculating multifunctionality from the set of functions indicated, for each "i" threshold 
# scaling can be by any function specified in "sc", i.e. max, mean, sd etc.
# here sc = maxx scales each function using the mean of the fives maximum values
# sc = qtle scales each function using the 0.95 quantile 

multi_all <- data.frame() 
for(i in seq(0.5, 0.8, by = 0.05)){
  multi <- multidiv(All_BigData[c("Aboveground_biomass", "Herbivory", "Pathogens", 
                                  "Plant_N_uptake", "Plant_P_uptake",
                                  "Soil_respiration", "Belowground_biomass", 
                                  "BGlucosidase", "Phosphatase", 
                                  "Carbon_storage"
  )]
  , threshold = i, sc = maxx
  )
  multi <- cbind(All_BigData["Plot_Nr"], multi)
  multi$Threshold <- rep(i,length(multi[,1]))
  multi_all = rbind(multi_all, multi)
}



MegaData <- merge(multi_all[c("Plot_Nr", "m", "Threshold")], 
                  All_BigData[c("Block", "Plot_Nr", "Species_richness", "Nitrogen", "Fungicide", "Combination", "CWM_SLA", "MPD_SLA_abundance", "MPD_SLA_presence")])

MegaData$Plot_Nr = as.factor(MegaData$Plot_Nr)
MegaData$Block = as.factor(MegaData$Block)
MegaData$Threshold_num = as.numeric(MegaData$Threshold)
MegaData$Threshold_num = scale(MegaData$Threshold_num)
# head(MegaData)




#### megamodel ####


# Main model

megamodel <- lmer(m ~ Threshold_num * (Nitrogen + Fungicide + Species_richness + CWM_SLA + MPD_SLA_presence)^2  
                  - CWM_SLA:MPD_SLA_presence - Threshold_num:CWM_SLA:MPD_SLA_presence 
                  + (1|Block) + (1|Combination) + 
                    (Threshold_num|Plot_Nr),
                  control = lmerControl(optimizer = "bobyqa"),
                  MegaData)

# model simplification, drops terms not significantly improving the overall model fit (likelihood-ratios).
stepba(megamodel)



