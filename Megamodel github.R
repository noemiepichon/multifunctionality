# 29/10/2021
# megamodel, code for the multifunctionality manuscript
# N.A. Pichon

library(lme4)

# load All_BigData
# load plot_description_scaled, which includes SLA and MPD values
# source the multidiv function
# source the stepba function

setwd("~/")
All_BigData <- read.table("All_BigData.txt", header = T)
plot_description_scaled <- read.table("plot_description_scaled.txt", header = T)
source("stepba function.R")
source("multidiv.R")


##### MegaData ####

# Calculating multifunctionality from the set of functions indicated, for each "i" threshold 
# scaling can be by any function specified in "sc", i.e. max, mean, sd etc.
# here sc = maxx scales each function using the mean of the fives maximum values
# sc = qtle scales each function using the 0.95 quantile 

multi_all <- data.frame() 
for(i in seq(0.5, 0.8, by = 0.05)){
  multi <- multidiv(All_BigData[c("Biomass", "CWM_Herb", "CWM_Path", 
                                  "N_ratio", "P_ratio",
                                  "Efflux", "Rootbiomass", 
                                  "BGlucosidase", "Phosphatase", 
                                  "carbon"
  )]
  , threshold = i, sc = maxx
  )
  multi <- cbind(All_BigData[1], multi)
  multi$Threshold <- rep(i,length(multi[,1]))
  multi_all = rbind(multi_all, multi)
}


MegaData <- merge(plot_description_scaled, multi_all[c("Plot_Nr", "m", "Threshold")])

MegaData$Plot_Nr = as.factor(MegaData$Plot_Nr)
MegaData$Block = as.factor(MegaData$Block)
MegaData$Threshold_num = as.numeric(MegaData$Threshold)
MegaData$Threshold_num = scale(MegaData$Threshold_num)
# head(MegaData)




#### megamodel ####


# Main model

megamodel <- lmer(m ~ Threshold_num * (Ni + Fz + SD + CS_SLA + MPD_SLA_pre)^2  
                  - CS_SLA:MPD_SLA_pre - Threshold_num:CS_SLA:MPD_SLA_pre 
                  + (1|Block) + (1|Comb) + 
                    (Threshold_num|Plot_Nr),
                  control = lmerControl(optimizer = "bobyqa"),
                  MegaData)

# model simplification, drops terms not significantly improving the overall model fit (likelihood-ratios).
stepba(megamodel)



