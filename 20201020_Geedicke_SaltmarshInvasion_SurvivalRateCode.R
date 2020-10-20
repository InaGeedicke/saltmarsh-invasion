##############################################################################
# Freshwater input drives invasion success of exotic plants in saltmarsh     #
#                                  communities                               #
#                                                                            #
#     Ina Geedicke, Anthony Manea, Jens Oldeland, Michelle R. Leishman       #
#                                                                            #
#                     Data analysis code, Austral Ecology                    #
##############################################################################


### SURVIVAL RATE ANALYSES -------------------------------------------------------

# Load packages -----------------------------------------------------------

library(dplyr)
library(survival) # survival functions (Kaplan Meyer Curve, Cox PH model)
library(survminer) # nice plotting with ggplot
library(gridExtra) # arranging of plots
library(multcomp) # pairwise comparison with Tukey/ Dunnet post hoc


# Load and clean data -----------------------------------------------------

sdata <- read.csv("20201020_Geedicke_SaltmarshInvasion_Survival_Exotics.csv", header = TRUE)
biodata <- read.csv("20201020_Geedicke_SaltmarshInvasion_Biomass.csv", header = TRUE)

# Summarize native saltmarsh biomass per pot (below and above ground and total) and add value to
# individual exotic plant per pot
NatAbove <- filter(biodata, Type == "native") %>% group_by(Pot) %>% dplyr::summarise(NatAB = sum(AboveBM))
NatBelow <- filter(biodata, Type == "native") %>% group_by(Pot) %>% dplyr::summarise(NatBB = sum(BelowBM))

newdat <- full_join(sdata, NatAbove, by = "Pot")
fsdata <- full_join(newdat, NatBelow, by = "Pot") %>% mutate(NatTotalB = NatBB + NatAB)

## Subset data per species
adata <- filter(fsdata, species == "Acetosa_saggitata")
bdata <- filter(fsdata, species == "Bidens_pilosa")
cdata <- filter(fsdata, species == "Conyza_parva")
pdata <- filter(fsdata, species == "Protasparagus_aethiopica")

head(fsdata)


# Kaplan-Meyer Survival Curve ---------------------------------------------

sp <- Surv(sdata$time, sdata$status) # calculate survival data: + stands for censored data
                                     # (alive) at the end of the trial

sp
sfit <- survfit(Surv(time, status)~ treatment + species, data=sdata) # Creates survival curve
summary(sfit) # gives summary of survival curve: survival in %, se & actual numbers of dead

fsfit <- survfit(Surv(time, status)~ treatment + species + NatTotalB, data=fsdata) # Creates survival curve with sm BM
summary(fsfit) # gives summary of survival curve: survival in %, se & actual numbers of dead

#survival curves per species
afit <- survfit(Surv(time, status)~ treatment + species, data=adata)
bfit <- survfit(Surv(time, status)~ treatment + species, data=bdata)
cfit <- survfit(Surv(time, status)~ treatment + species, data=cdata)
pfit <- survfit(Surv(time, status)~ treatment + species, data=pdata)

# general plotting using survminer package
ggsurvplot(pfit, conf.int=FALSE, pval=TRUE, 
                legend.labs=c("HSHN", "HSLN/control", "LSHN", "LSLN"), legend.title="Treatment",  
                #palette=c("black", "darkgrey", "black", "darkgrey" ), xlab = "time [weeks]",
                palette = c("turquoise4", "darkred", "turquoise4", "darkred" ),
                ylab = expression("survival probability "~italic("Protasparagus aethiopica")), 
                linetype = c('solid', 'solid', 'dashed', 'dashed'))


# Plot survival curves for each species in a grid
splots <- list()
splots[[1]] <- ggsurvplot(afit, data = adata, conf.int = FALSE, pval = FALSE, legend = "none",
                          #palette = c("turquoise", "darkred", "turquoise", "darkred" ),
                          palette = c("turquoise4", "darkred", "turquoise4", "darkred" ),
                          xlab = NULL, 
                          ylab = "survival probability ", 
                          linetype = c('solid', 'solid', 'dashed', 'dashed'))

splots[[2]] <- ggsurvplot(bfit, data = bdata, pval = FALSE, legend = "none", 
                          #palette = c("black", "darkgrey", "black", "darkgrey" ),
                          palette = c("turquoise4", "darkred", "turquoise4", "darkred" ),
                          xlab = "time [weeks]", 
                          ylab = "survival probability ", 
                          linetype = c('solid', 'solid', 'dashed', 'dashed'))

splots[[3]] <- ggsurvplot(cfit, data = cdata, pval = FALSE, legend = "none", 
                          #palette = c("black", "darkgrey", "black", "darkgrey" ), 
                          palette = c("turquoise4", "darkred", "turquoise4", "darkred" ),
                          xlab = NULL, 
                          ylab = NULL, 
                          linetype = c('solid', 'solid', 'dashed', 'dashed'))

splots[[4]] <- ggsurvplot(pfit, data = pdata, pval = FALSE, legend = "none", 
                          #palette = c("black", "darkgrey", "black", "darkgrey" ), 
                          palette = c("turquoise4", "darkred", "turquoise4", "darkred" ),
                          xlab = "time [weeks]", 
                          ylab = NULL, 
                          linetype = c('solid', 'solid', 'dashed', 'dashed'))

arrange_ggsurvplots(splots,
                    ncol = 2, nrow = 2)


# Cox PH Model ------------------------------------------------------------

# Cox PH regression can assess the effect of both categorical and continuous variables, and can 
# model the effect of multiple variables at once.
# Cox PH regression models the natural log of the hazard at time t, denoted h(t), as a 
# function of the baseline hazard (h0(t)) (the hazard for an individual where all exposure
# variables are 0) and multiple exposure variables x1, x2,..., xp. The form of the Cox PH model 
# is:

#  log(h(t))=log(h0(t))+β1x1+β2x2+...+βpxp

# If you exponentiate both sides of the equation, and limit the right hand side to just a single
# categorical exposure variable (x1) with two groups (x1=1 for exposed and x1=0 for unexposed), 
# the equation becomes:
  
#  h1(t)=h0(t)+eβ1x1

# Rearranging that equation lets you estimate the hazard ratio, comparing the exposed to the
#  unexposed individuals at time t:
  
#  HR(t)= h0(t)eβ1 / h0(t) = eβ1

# This model shows that the hazard ratio is eβ1, and remains constant over time t (hence 
# the name proportional hazards regression). The β values are the regression coefficients that
# are estimated from the model, and represent the log(HazardRatio) for each unit increase in the
# corresponding predictor variable. The interpretation of the hazards ratio depends on the 
# measurement scale of the predictor variable, but in simple terms, a positive coefficient indicates
# worse survival and a negative coefficient indicates better survival for the variable in question.

# exp(coef) contains eβ1, which is the hazard ratio (HR)
# HR=1: No effect
# HR>1: Increase in hazard
# HR<1: Reduction in hazard (protective)

# Source: http://bioconnector.org/workshops/r-survival.html#cox_ph_model

fit <- coxph(Surv(time, status)~ treatment + species, data=sdata)
fit


# Multiple comparison -----------------------------------------------------

# General Linear Hypothesis (glht)
# Dunnett: multiple comparison if you have a control group, here HSLN
# Tukey: multiple comparison without control


## Multiple comparison Total

# control group must be first > rearrange data, so that HSLN is first group in col treatment
library(gdata)
target <- c("HSLN", "HSHN", "LSHN", "LSLN")
# for sdata
sdata$treatment <- reorder.factor(sdata$treatment, new.order=target)
sdatasort <- sdata %>% arrange(treatment)

# Again for fsdata
fsdata$treatment <- reorder.factor(fsdata$treatment, new.order=target)
fsdatasort <- fsdata %>% arrange(treatment)


# Multiple comparison for the total dataset comparing species and treatment
fitsort <- coxph(Surv(time, status)~ treatment + species, data=sdatasort) # Cox PH Model
glht_fit <- glht(model = fitsort, linfct = mcp(treatment = "Dunnett", species = "Tukey"))
summary(glht_fit)

# Multiple comparison for the total dataset comparing species and treatment, including native saltmarsh biomass
ffitsort <- coxph(Surv(time, status)~ treatment + species + NatTotalB, data=fsdatasort) # Cox PH Model
glht_fit <- glht(model = ffitsort, linfct = mcp(treatment = "Dunnett", species = "Tukey", NatTotalB = "Tukey"))
summary(glht_fit)


## Multiple comparison per exotic species

# Filter sorted dataset per species
asdatasort <- filter(sdatasort, species == "Acetosa_saggitata")
bsdatasort <- filter(sdatasort, species == "Bidens_pilosa")
csdatasort <- filter(sdatasort, species == "Conyza_parva")
psdatasort <- filter(sdatasort, species == "Protasparagus_aethiopica")

# Multiple comparison per species against treatment
afitsort <- coxph(Surv(time, status)~ treatment, data=asdatasort) # Cox PH Model
aglht_fit <- glht(model = afitsort, linfct = mcp(treatment = "Dunnett"))
summary(aglht_fit)

bfitsort <- coxph(Surv(time, status)~ treatment, data=bsdatasort) # Cox PH Model
bglht_fit <- glht(model = bfitsort, linfct = mcp(treatment = "Dunnett"))
summary(bglht_fit)

cfitsort <- coxph(Surv(time, status)~ treatment, data=csdatasort) # Cox PH Model
cglht_fit <- glht(model = cfitsort, linfct = mcp(treatment = "Dunnett"))
summary(cglht_fit)

pfitsort <- coxph(Surv(time, status)~ treatment, data=psdatasort) # Cox PH Model
pglht_fit <- glht(model = pfitsort, linfct = mcp(treatment = "Dunnett"))
summary(pglht_fit)


# Multiple comparison per species against treatment using TUKEY
afitsort <- coxph(Surv(time, status)~ treatment, data=asdatasort) # Cox PH Model
aglht_fit <- glht(model = afitsort, linfct = mcp(treatment = "Tukey"))
summary(aglht_fit)

bfitsort <- coxph(Surv(time, status)~ treatment, data=bsdatasort) # Cox PH Model
bglht_fit <- glht(model = bfitsort, linfct = mcp(treatment = "Tukey"))
summary(bglht_fit)

cfitsort <- coxph(Surv(time, status)~ treatment, data=csdatasort) # Cox PH Model
cglht_fit <- glht(model = cfitsort, linfct = mcp(treatment = "Tukey"))
summary(cglht_fit)

pfitsort <- coxph(Surv(time, status)~ treatment, data=psdatasort) # Cox PH Model
pglht_fit <- glht(model = pfitsort, linfct = mcp(treatment = "Tukey"))
summary(pglht_fit)
