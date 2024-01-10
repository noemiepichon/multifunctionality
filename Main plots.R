# 10/01/2024
# how to plot the main figure, code for Yao Xiao and others
# N.A. Pichon

library(wesanderson)
library(ggplot2)
library(cowplot)


source("Megamodel github.R") # or run the code to get the MegaData dataset, don't need to run stepba


Final_model <- lmer(formula = m ~ Threshold_num + Nitrogen + Fungicide + Species_richness + 
                      CWM_SLA + MPD_SLA_presence + (1 | Block) + (1 | Combination) + 
                      (Threshold_num | Plot_Nr) + Nitrogen:MPD_SLA_presence + Fungicide:MPD_SLA_presence + 
                      Species_richness:CWM_SLA, data = MegaData)





# Model prediction Species_richness x SLA quantiles

newdat <- expand.grid(
  Species_richness=c(-1.05734693249371,  -0.348699520290478, 0.596163695980495 , 3.43075334479342),
  # CWM_SLA=quantile(MegaData$CWM_SLA, names = F),
  CWM_SLA = c(-1.421408, -0.7263824, -0.1030613,  0.5692412, 1.40453), # only middle quantiles, and two more, for 8 species test
  Nitrogen = 0, Fungicide = 0, MPD_SLA_presence = 0,
  Threshold_num = 0,
  m = 0
)

newdat = newdat[-c(4,20),]

# calculate errors
newdat$m <- predict(Final_model, newdat, re.form=NA)
mm <- model.matrix(terms(Final_model),newdat)
## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(Final_model), mm))
tvar1 <- pvar1+VarCorr(Final_model)$Block[1]+VarCorr(Final_model)$Block[1]+
  VarCorr(Final_model)$Plot_Nr[1]  ## must be adapted for more complex models
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , plo = newdat$m-cmult*sqrt(pvar1)
  , phi = newdat$m+cmult*sqrt(pvar1)
  , tlo = newdat$m-cmult*sqrt(tvar1)
  , thi = newdat$m+cmult*sqrt(tvar1)
)

# backtransform
newdat["Species_richness"] <- pd_attributes_variable$`scaled:center`["Species_richness"] + pd_attributes_variable$`scaled:scale`["Species_richness"]*newdat["Species_richness"]
newdat["CWM_SLA"] <- round(pd_attributes_variable$`scaled:center`["CWM_SLA"] + pd_attributes_variable$`scaled:scale`["CWM_SLA"]*newdat["CWM_SLA"],1)

gdbp <- wes_palette("GrandBudapest1", 7, type = "continuous")

# Model Species_richness x SLA multi / PM - Species_richness x SLA x Thr no points
ggplot(newdat, aes(Species_richness, m, col=as.factor(CWM_SLA))) + 
  # facet_grid(~Threshold_num) +
  geom_ribbon(aes(ymin = plo, ymax = phi), linetype = 0, alpha = 0.1)+
  geom_smooth(method = lm, se = F, size = 1.5) +
  theme_cowplot()+
  coord_cartesian(ylim = c(0.3, 1))+
  coord_cartesian(ylim = c(0.4, 0.75))+
  scale_x_continuous(breaks=c(1,4,8,20))+
  labs(x = "Species richness", y = "Multifunctionality")+
  scale_colour_manual(values = gdbp, name = "SLA quantiles", guide = guide_legend(reverse = TRUE))

#

