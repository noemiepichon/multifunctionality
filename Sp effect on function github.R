# 04/12/2023 - updated EcoLetters
# sp effect on function
# code for multifunctionality
# N.A. Pichon


library(tidyr)

# setwd("~/")
Mean_SLA <- read.table("Mean_SLA.txt", header = T)
All_BigData <- read.table("All_BigData.txt", header = T)
sp_presence <- read.table("sp_presence.txt", header = T)
source("beta_div_functions.R")
#


#### Prepare dataset ####


All_BigData$Aboveground_biomass = sqrt(All_BigData$Aboveground_biomass)
All_BigData$Herbivory = sqrt(All_BigData$Herbivory)
# All_BigData$Pathogens = All_BigData$Pathogens # no transfomration needed
All_BigData$Soil_respiration = log(All_BigData$Soil_respiration)
All_BigData$Plant_N_uptake = log(All_BigData$Plant_N_uptake)
All_BigData$Plant_P_uptake = log(All_BigData$Plant_P_uptake)
All_BigData$Belowground_biomass = log(All_BigData$Belowground_biomass)
All_BigData$BGlucosidase = sqrt(All_BigData$BGlucosidase)
All_BigData$Phosphatase = log(All_BigData$Phosphatase)
All_BigData$Carbon_storage = sqrt(All_BigData$Carbon_storage+1)


All_BigData_presence = merge(All_BigData, sp_presence)
All_BigData_presence[,c("Species_richness", "Nitrogen", "Fungicide", "CWM_SLA", "MPD_SLA_abundance", "MPD_SLA_presence", 
                        "Aboveground_biomass", "Herbivory", "Pathogens", "Soil_respiration", "Plant_N_uptake", "Plant_P_uptake", "Belowground_biomass", "BGlucosidase", "Phosphatase", 
                        "Carbon_storage")] = scale(All_BigData_presence[,c("Species_richness", "Nitrogen", "Fungicide", "CWM_SLA", "MPD_SLA_abundance", "MPD_SLA_presence", 
                                                                           "Aboveground_biomass", "Herbivory", "Pathogens", "Soil_respiration","Plant_N_uptake", 
                                                                           "Plant_P_uptake", "Belowground_biomass", "BGlucosidase", "Phosphatase", "Carbon_storage")])



# Extract the effect of each species for each function using lm

Sp_Effect_pre = data.frame()
for(i in names(All_BigData_presence[,c(11:20)])){
  mod = lm(All_BigData_presence[,i] ~ Nitrogen*Fungicide + Am + Ao + As + Be + Cb + Cj + Dc + Dg + Fr + Ga + 
             Hl + Hp + Hs + Lp + Pg + Pm + Pt + Ra + Sp + To + 
             Block, All_BigData_presence)
  Sp_df <- as.data.frame(summary(mod)$coefficients[,c(1,2,4)])
  Sp_df$Function <- rep(i,length(Sp_df[,1]))
  Variable <- rownames(Sp_df)
  rownames(Sp_df) <- NULL
  Sp_df <- cbind(Variable,Sp_df)
  Sp_Effect_pre = rbind(Sp_Effect_pre, Sp_df)
}

## Here is the change
# Sp_Effect_pre$Estimate[Sp_Effect_pre$`Pr(>|t|)` >0.05] <- 0
# Sp_Effect_pre$`Std. Error`[Sp_Effect_pre$`Pr(>|t|)` >0.05] <- 0
Sp_Effect_pre$Estimate = Sp_Effect_pre$Estimate/Sp_Effect_pre$`Std. Error`


Sp_Effect_pre <- Sp_Effect_pre[!Sp_Effect_pre$Variable=="(Intercept)",]
Sp_Effect_pre <- Sp_Effect_pre[!Sp_Effect_pre$Variable=="Nitrogen:Fungicide",]
Sp_Effect_pre <- Sp_Effect_pre[!Sp_Effect_pre$Variable=="Fungicide",]
Sp_Effect_pre <- Sp_Effect_pre[!Sp_Effect_pre$Variable=="Nitrogen",]
Sp_Effect_pre <- Sp_Effect_pre[!Sp_Effect_pre$Variable=="Block",]





#### Extract the effect of mean and variance in SLA ####



RowVar <-  function(x, ...){
  rowSums((x-rowMeans(x, ...))^2, ...)/(dim(x)[2]-1)
}



##
## positive
##

chao_pre_pos <- subset(Sp_Effect_pre, Estimate >0)
chao_pre_pos_test = spread(chao_pre_pos[,c("Variable", "Function", "Estimate")], Function, Estimate)
rownames(chao_pre_pos_test) = (chao_pre_pos_test$Variable)
chao_pre_pos_test$Variable <- NULL
chao_pre_pos_test[is.na(chao_pre_pos_test)] <- 0


chao_pre_pos_scol = t(chao_pre_pos_test)
sp_pairs = combn(ncol(chao_pre_pos_scol), 4)

df_pre_pos_scol = data.frame()
for(i in c(1:ncol(sp_pairs))){
  yep = CqN(chao_pre_pos_scol[,(sp_pairs[,i])], from = 0, to = 2, interval = 1)[,c("q", "qD.alpha", "CqN")]
  yep$Sp1 = colnames(chao_pre_pos_scol[,(sp_pairs[,i])])[1]
  yep$Sp2 = colnames(chao_pre_pos_scol[,(sp_pairs[,i])])[2]
  yep$Sp3 = colnames(chao_pre_pos_scol[,(sp_pairs[,i])])[3]
  yep$Sp4 = colnames(chao_pre_pos_scol[,(sp_pairs[,i])])[4]
  df_pre_pos_scol = rbind(df_pre_pos_scol, yep)
}

df_pre_pos_scol = merge(df_pre_pos_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp1", by.y = "Species")
df_pre_pos_scol = merge(df_pre_pos_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp2", by.y = "Species")
df_pre_pos_scol = merge(df_pre_pos_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp3", by.y = "Species")
df_pre_pos_scol = merge(df_pre_pos_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp4", by.y = "Species")
head(df_pre_pos_scol)

df_pre_pos_scol$Mean_SLA = rowMeans(df_pre_pos_scol[,c(8:11)])
df_pre_pos_scol$Var_SLA = RowVar(df_pre_pos_scol[,c(8:11)])



# # visualise the output
# plot(CqN ~ Mean_SLA, subset(df_pre_pos_scol, q == 2))
# abline(lm(CqN ~ Mean_SLA, subset(df_pre_pos_scol, q == 2)))



#


#



##
## negative
##

chao_pre_neg <- subset(Sp_Effect_pre, Estimate <0)
chao_pre_neg_test = spread(chao_pre_neg[,c("Variable", "Function", "Estimate")], Function, Estimate)
rownames(chao_pre_neg_test) = (chao_pre_neg_test$Variable)
chao_pre_neg_test$Variable <- NULL
chao_pre_neg_test[is.na(chao_pre_neg_test)] <- 0


chao_pre_neg_scol = t(chao_pre_neg_test)
sp_pairs = combn(ncol(chao_pre_neg_scol), 4)

df_pre_neg_scol = data.frame()
for(i in c(1:ncol(sp_pairs))){
  yep = CqN(chao_pre_neg_scol[,(sp_pairs[,i])], from = 0, to = 2, interval = 1)[,c("q", "qD.alpha", "CqN")]
  yep$Sp1 = colnames(chao_pre_neg_scol[,(sp_pairs[,i])])[1]
  yep$Sp2 = colnames(chao_pre_neg_scol[,(sp_pairs[,i])])[2]
  yep$Sp3 = colnames(chao_pre_neg_scol[,(sp_pairs[,i])])[3]
  yep$Sp4 = colnames(chao_pre_neg_scol[,(sp_pairs[,i])])[4]
  df_pre_neg_scol = rbind(df_pre_neg_scol, yep)
}

df_pre_neg_scol = merge(df_pre_neg_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp1", by.y = "Species")
df_pre_neg_scol = merge(df_pre_neg_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp2", by.y = "Species")
df_pre_neg_scol = merge(df_pre_neg_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp3", by.y = "Species")
df_pre_neg_scol = merge(df_pre_neg_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp4", by.y = "Species")
head(df_pre_neg_scol)


df_pre_neg_scol$Mean_SLA = rowMeans(df_pre_neg_scol[,c(8:11)])
df_pre_neg_scol$Var_SLA = RowVar(df_pre_neg_scol[,c(8:11)])



# # visualise the output
# plot(CqN ~ Mean_SLA, subset(df_pre_neg_scol, q == 2))
# abline(lm(CqN ~ Mean_SLA, subset(df_pre_neg_scol, q == 2)))



#




#### Test the effect of mean and variance in SLA on functional effects, by randomisation ####


# Using only q == 2, which gives more weight to the large functional effects




# positive

df_pre_pos_scol2 <- subset(df_pre_pos_scol, q == 2)[,c(7,12,13)]
head(df_pre_pos_scol2)


Pre_pos_scol_random = data.frame()
for(i in c(1:1000)){
  df_pre_pos_scol2[,c(2:3)] = df_pre_pos_scol2[sample(nrow(df_pre_pos_scol2)), c(2:3)]
  mod = lm(CqN ~ Mean_SLA + Var_SLA, df_pre_pos_scol2)
  yep <- as.data.frame(summary(mod)$coefficients[,c(1,2)])
  Variprele <- rownames(yep)
  rownames(yep) <- NULL
  yep <- cbind(Variprele,yep)
  Pre_pos_scol_random = rbind(Pre_pos_scol_random, yep)
}

Pre_pos_scol_random_meanSLA = subset(Pre_pos_scol_random, Variprele == "Mean_SLA")
# mean(Pre_pos_scol_random_meanSLA$Estimate)
# sd(Pre_pos_scol_random_meanSLA$Estimate)
quantile(Pre_pos_scol_random_meanSLA$Estimate, c(0.025, 0.975))



# visualise the output
perm_pre_pos_scol <- lm(CqN ~ Mean_SLA + Var_SLA, subset(df_pre_pos_scol, q == 2))

hist(Pre_pos_scol_random_meanSLA$Estimate, xlim = c(-0.01,0.01), main = "Positive, mean SLA", xlab = "Estimate")
abline(v = quantile(Pre_pos_scol_random_meanSLA$Estimate, 0.025))
abline(v = quantile(Pre_pos_scol_random_meanSLA$Estimate, 0.975))
abline(v = as.data.frame(perm_pre_pos_scol$coefficients['Mean_SLA'])[,1], col = "red")

Pre_pos_scol_random_varSLA = subset(Pre_pos_scol_random, Variprele == "Var_SLA")
hist(Pre_pos_scol_random_varSLA$Estimate, xlim = c(-0.0005,0.0005), main = "Positive, var SLA", xlab = "Estimate")
abline(v = quantile(Pre_pos_scol_random_varSLA$Estimate, 0.025))
abline(v = quantile(Pre_pos_scol_random_varSLA$Estimate, 0.975))
abline(v = as.data.frame(perm_pre_pos_scol$coefficients['Var_SLA'])[,1], col = "red")

#

####
# Positive values, effect of SLA mean and variance
perm_pre_pos_scol$coefficients
# quantiles mean 
quantile(Pre_pos_scol_random_meanSLA$Estimate, c(0.025, 0.975))
# quantiles variance
quantile(Pre_pos_scol_random_varSLA$Estimate, c(0.025, 0.975))
####



# negative

df_pre_neg_scol2 <- subset(df_pre_neg_scol, q == 2)[,c(7,12,13)]
head(df_pre_neg_scol2)


Pre_neg_scol_random = data.frame()
for(i in c(1:1000)){
  df_pre_neg_scol2[,c(2:3)] = df_pre_neg_scol2[sample(nrow(df_pre_neg_scol2)), c(2:3)]
  mod = lm(CqN ~ Mean_SLA + Var_SLA, df_pre_neg_scol2)
  yep <- as.data.frame(summary(mod)$coefficients[,c(1,2)])
  Variprele <- rownames(yep)
  rownames(yep) <- NULL
  yep <- cbind(Variprele,yep)
  Pre_neg_scol_random = rbind(Pre_neg_scol_random, yep)
}


Pre_neg_scol_random_meanSLA = subset(Pre_neg_scol_random, Variprele == "Mean_SLA")
# mean(Pre_neg_scol_random_varSLA$Estimate)
# sd(Pre_neg_scol_random_varSLA$Estimate)
# quantile(Pre_neg_scol_random_varSLA$Estimate, c(0.025, 0.975))



# visualise the output
perm_pre_neg_scol <- lm(CqN ~ Mean_SLA + Var_SLA, subset(df_pre_neg_scol, q == 2))


hist(Pre_neg_scol_random_meanSLA$Estimate, xlim = c(-0.005,0.005), main = "Negative, mean SLA", xlab = "Estimate")
abline(v = quantile(Pre_neg_scol_random_meanSLA$Estimate, 0.025))
abline(v = quantile(Pre_neg_scol_random_meanSLA$Estimate, 0.975))
abline(v = as.data.frame(perm_pre_neg_scol$coefficients['Mean_SLA'])[,1], col = "red")

Pre_neg_scol_random_varSLA = subset(Pre_neg_scol_random, Variprele == "Var_SLA")
hist(Pre_neg_scol_random_varSLA$Estimate, xlim = c(-0.001,0.001), main = "Negative, var SLA", xlab = "Estimate")
abline(v = quantile(Pre_neg_scol_random_varSLA$Estimate, 0.025))
abline(v = quantile(Pre_neg_scol_random_varSLA$Estimate, 0.975))
abline(v = as.data.frame(perm_pre_neg_scol$coefficients['Var_SLA'])[,1], col = "red")


####
# Negative values, effect of SLA mean and variance
perm_pre_neg_scol$coefficients
# quantiles mean 
quantile(Pre_neg_scol_random_meanSLA$Estimate, c(0.025, 0.975))
# quantiles variance
quantile(Pre_neg_scol_random_varSLA$Estimate, c(0.025, 0.975))
####





#### Test the effect of mean and variance in SLA on trade-offs, by randomisation ####





chao_pre_to <- subset(Sp_Effect_pre, !Estimate == 0)
chao_pre_to_test = spread(chao_pre_to[,c("Variable", "Function", "Estimate")], Function, Estimate)
rownames(chao_pre_to_test) = (chao_pre_to_test$Variable)
chao_pre_to_test$Variable <- NULL
chao_pre_to_test[is.na(chao_pre_to_test)] <- 0

#
# 
#

chao_pre_to_scol = t(chao_pre_to_test)
sp_pairs = combn(ncol(chao_pre_to_scol), 2)

df_pre_to_scol = data.frame()
for(i in c(1:ncol(sp_pairs))){
  yep = as.data.frame(chao_pre_to_scol[,(sp_pairs[,i])])
  yep$prod = yep[,1]*yep[,2]
  yep$to <- ifelse(yep$prod <0, "1", ifelse(Sp_Effect_pre$Estimate >0, "0", "0"))
  yep2 = as.data.frame(sum(as.numeric(yep$to)))
  yep2$Sp1 = colnames(chao_pre_to_scol[,(sp_pairs[,i])])[1]
  yep2$Sp2 = colnames(chao_pre_to_scol[,(sp_pairs[,i])])[2]
  df_pre_to_scol = rbind(df_pre_to_scol, yep2)
}


names(df_pre_to_scol) = c("Trade_offs", "Sp1", "Sp2")

df_pre_to_scol = merge(df_pre_to_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp1", by.y = "Species")
df_pre_to_scol = merge(df_pre_to_scol, Mean_SLA[,c("CWM_SLA", "Species")], by.x = "Sp2", by.y = "Species")
head(df_pre_to_scol)

df_pre_to_scol$Mean_SLA = rowMeans(df_pre_to_scol[,c(4:5)])
df_pre_to_scol$Var_SLA = RowVar(df_pre_to_scol[,c(4:5)])


mod_pre_to_scol <- lm(Trade_offs ~ Mean_SLA + Var_SLA, df_pre_to_scol)
mod_pre_to_scol$coefficients
#





df_pre_to_scol2 <- df_pre_to_scol[,c(3,6,7)]
head(df_pre_to_scol2)


Pre_pos_to_random = data.frame()
for(i in c(1:1000)){
  df_pre_to_scol2[,c(2:3)] = df_pre_to_scol2[sample(nrow(df_pre_to_scol2)), c(2:3)]
  mod = lm(Trade_offs ~ Mean_SLA + Var_SLA, df_pre_to_scol2)
  yep <- as.data.frame(summary(mod)$coefficients[,c(1,2)])
  Variable <- rownames(yep)
  rownames(yep) <- NULL
  yep <- cbind(Variable,yep)
  Pre_pos_to_random = rbind(Pre_pos_to_random, yep)
}

head(Pre_pos_to_random)

Pre_pos_to_random_meanSLA = subset(Pre_pos_to_random, Variable == "Mean_SLA")
hist(Pre_pos_to_random_meanSLA$Estimate, main = "Trade-offs, mean SLA", xlab = "Estimate")#, xlim = c(-0.02,0.02)
abline(v = quantile(Pre_pos_to_random_meanSLA$Estimate, 0.025))
abline(v = quantile(Pre_pos_to_random_meanSLA$Estimate, 0.975))
abline(v = as.data.frame(mod_pre_to_scol$coefficients['Mean_SLA'])[,1], col = "red")
quantile(Pre_pos_to_random_meanSLA$Estimate, c(0.025, 0.975))


Pre_pos_to_random_varSLA = subset(Pre_pos_to_random, Variable == "Var_SLA")
hist(Pre_pos_to_random_varSLA$Estimate, main = "Trade-offs, var SLA", xlab = "Estimate")#, xlim = c(-0.02,0.02))
abline(v = quantile(Pre_pos_to_random_varSLA$Estimate, 0.025))
abline(v = quantile(Pre_pos_to_random_varSLA$Estimate, 0.975))
abline(v = as.data.frame(mod_pre_to_scol$coefficients['Var_SLA'])[,1], col = "red")
quantile(Pre_pos_to_random_varSLA$Estimate, c(0.025, 0.975))
#
#




####
# Trade-offs, effect of SLA mean and variance
mod_pre_to_scol$coefficients
# quantiles mean 
quantile(Pre_pos_to_random_meanSLA$Estimate, c(0.025, 0.975))
# quantiles variance
quantile(Pre_pos_to_random_varSLA$Estimate, c(0.025, 0.975))
####
