# 10/01/2024
# how to plot the main figure, code for Yao Xiao and others
# N.A. Pichon

library(wesanderson)
library(ggplot2)
library(cowplot)


# source("Megamodel github.R") # or run the code to get the MegaData dataset, don't need to run stepba


Final_model <- lmer(formula = m ~ Threshold_num + Nitrogen + Fungicide + Species_richness + 
                      CWM_SLA + MPD_SLA_presence + (1 | Block) + (1 | Combination) + 
                      (Threshold_num | Plot_Nr) + Nitrogen:MPD_SLA_presence + Fungicide:MPD_SLA_presence + 
                      Species_richness:CWM_SLA, control = lmerControl(optimizer = "bobyqa"), data = MegaData)






## interaction figures

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





# Model prediction MPD_SLA_presence x Nitrogen

newdat <- expand.grid(
  MPD_SLA_presence = MegaData$MPD_SLA_presence, 
  Species_richness = 0,
  Nitrogen=c(-0.9985, 0.9985), 
  CWM_SLA = 0, Fungicide = 0,
  Threshold_num = 0,
  m = 0
)

newdat$m <- predict(Final_model, newdat, re.form=NA)
mm <- model.matrix(terms(Final_model),newdat)
## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(Final_model),mm))
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
newdat["MPD_SLA_presence"] <- (pd_attributes_variable$`scaled:center`["MPD_SLA_presence"] + pd_attributes_variable$`scaled:scale`["MPD_SLA_presence"]*newdat["MPD_SLA_presence"])

# Model MPD_SLA_ab x Nitrogen multi
ggplot(newdat, aes(MPD_SLA_presence, m, col=as.factor(Nitrogen))) + 
  # facet_grid(~Threshold_num) +
  geom_ribbon(aes(ymin = plo, ymax = phi), linetype = 0, alpha = 0.1)+
  geom_line(aes(y=phi), size = 0.1) +
  geom_smooth(method = lm, se = F, size = 1.5) +
  geom_line(aes(y=plo), size = 0.1) +
  theme_cowplot()+
  coord_cartesian(ylim = c(0.4, 0.75))+
  labs(x = "Functional diversity (MPD)", y = "Multifunctionality")+
  scale_colour_manual(values = c("#99CCFF", "#336699"), name = "Treatment", labels = c("without nitrogen", "with nitrogen"))
#




# Model prediction MPD_SLA_presence x Fungicide 

newdat <- expand.grid(
  MPD_SLA_presence = MegaData$MPD_SLA_presence,
  Species_richness = 0,
  Nitrogen=0, 
  CWM_SLA = 0, Fungicide = c(-0.9985, 0.9985),
  Threshold_num = 0,#c(-1.4995039, 0, 1.4995039),#seq(0.5,0.8,0.05),#0.65,
  m = 0
)
newdat$m <- predict(Final_model, newdat, re.form=NA)
mm <- model.matrix(terms(Final_model),newdat)
## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(Final_model),mm))
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
newdat["MPD_SLA_presence"] <- (pd_attributes_variable$`scaled:center`["MPD_SLA_presence"] + pd_attributes_variable$`scaled:scale`["MPD_SLA_presence"]*newdat["MPD_SLA_presence"])

# Model MPD_SLA_presence x Nitrogen
ggplot(newdat, aes(MPD_SLA_presence, m, col=as.factor(Fungicide))) + 
  # facet_grid(~Threshold_num) +
  geom_ribbon(aes(ymin = plo, ymax = phi), linetype = 0, alpha = 0.1)+
  geom_line(aes(y=plo), size = 0.1) +
  geom_line(aes(y=phi), size = 0.1) +
  geom_smooth(method = lm, se = F, size = 1.5) +
  theme_cowplot()+
  coord_cartesian(ylim = c(0.4, 0.75))+
  labs(x = "Functional diversity (MPD)", y = "Multifunctionality")+
  scale_colour_manual(values = c("#FF6666","#990000"), name = "Treatment", labels = c("control", "with fungicide"))
#











## Main figure


summary(Final_model)
MegaModelSummary <- as.data.frame(coef(summary(Final_model)))
MegaModelSummary$Effect <- rownames(MegaModelSummary)

# confint function
conf_intervals <- confint(Final_model, method="boot", level = 0.95) # this eats time
conf_intervals2 = as.data.frame(conf_intervals[c(7:16),])
names(conf_intervals2) = c("intlow", "inthigh")
conf_intervals2$Effect <- rownames(conf_intervals2)
MegaModelSummary <- merge(MegaModelSummary, conf_intervals2)

MegaModelSummary$Order <- c(rep(1,3), 2,1,1,2,1,2,1)#c(rep(1,6), rep(2,3))

MegaModelSummary$Effect <- factor(MegaModelSummary$Effect, levels = MegaModelSummary$Effect[order(MegaModelSummary$Order,
                                                                                                  abs(MegaModelSummary$Estimate),
                                                                                                  decreasing = c(T, F),
                                                                                                  method="radix")])

ggplot(subset(MegaModelSummary, Estimate<0.1 & Estimate>-0.1), aes(x = Effect, y = Estimate, fill = (Order)))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin=intlow, ymax=inthigh, color = "black"), width=.3, position = position_dodge(width=0.5), size = 0.4)+ # confint
  geom_hline(yintercept = 0)+
  theme_cowplot()+
  theme(legend.position = "")+
  coord_flip()+
  scale_x_discrete(labels=c("Fungicide x MPD", "Nitrogen x MPD", "Species richness x SLA",
                            "Functional diversity (MPD)", "Fungicide", "Species richness", "Nitrogen", "Specific leaf area (SLA)"), name = c())
#
