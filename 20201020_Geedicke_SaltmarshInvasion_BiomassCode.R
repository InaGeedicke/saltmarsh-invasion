##############################################################################
# Freshwater input drives invasion success of exotic plants in saltmarsh     #
#                                  communities                               #
#                                                                            #
#     Ina Geedicke, Anthony Manea, Jens Oldeland, Michelle R. Leishman       #
#                                                                            #
#                     Data analysis code, Austral Ecology                    #
##############################################################################


### BIOMASS ANALYSES -------------------------------------------------------


# Load packages -----------------------------------------------------------

library(plyr)
library(dplyr)     # cleaning of data
library(stargazer) # creates tables of model results
library(ggplot2)   # plotting
library(car)       # Anova()
library(multcomp)  # post hoc test (glht > Dunnett)
library(pscl)      # pseudo R square for glm
library(effects)   # explore effects of GLM


# Load and clean data -----------------------------------------------------

biodata <- read.csv("20201020_Geedicke_SaltmarshInvasion_Biomass.csv", header = TRUE)

# Summarize native saltmarsh biomass per pot (below and above ground and total) and add value to
# individual exotic plant per pot
NatAbove <- filter(biodata, Type == "native") %>% group_by(Pot) %>% dplyr::summarise(NatAB = sum(AboveBM))
NatBelow <- filter(biodata, Type == "native") %>% group_by(Pot) %>% dplyr::summarise(NatBB = sum(BelowBM))

newdat <- full_join(biodata, NatAbove, by = "Pot")
finaldat <- full_join(newdat, NatBelow, by = "Pot") %>% mutate(NatTotalB = NatBB + NatAB) %>% 
                                                        mutate(ExTotalB = AboveBM + BelowBM)  %>%
                                                        filter(Type == "exotic")

## Summarize data to include native saltmarsh only

SMdat <- filter(biodata, Type == "native") %>% mutate(NatTotalB = AboveBM + BelowBM)


#  Generalized Linear Model (GLM) to compare biomass and treatments -------

## create factorial column for native saltmarsh biomass (high and low BM)
finaldat$fNatB <- as.factor(ifelse(finaldat$NatTotalB > 10000, "High", "Low"))

## All native species
glm.SMall <- glm(NatTotalB ~ Treatment + Species_short, family="Gamma"(link="log"), 
                 data=SMdat)
summary(glm.SMall)
Anova(glm.SMall)
Tukey.testSMall <- glht(glm.SMall, linfct=mcp(Treatment="Tukey")) 
summary(Tukey.testSMall)

## All exotic species
glm.EXall <- glm(ExTotalB ~ Treatment + fNatB + Species_short, family="Gamma"(link="log"), 
               data=finaldat)
summary(glm.EXall)
Anova(glm.EXall)
Tukey.testEXall <- glht(glm.EXall, linfct=mcp(Treatment="Tukey")) 
summary(Tukey.testEXall)

## Protasparagus aethiopicus
glm.PA<-glm(ExTotalB ~ Treatment + fNatB, family="Gamma"(link="log"), 
            data=finaldat, subset=Species_short=="PA")
summary(glm.PA)
Anova(glm.PA)
Tukey.testPA <- glht(glm.PA, linfct=mcp(Treatment="Tukey")) 
summary(Tukey.testPA)

## Acetosa sagittata
glm.AS<-glm(ExTotalB ~ Treatment + fNatB, family="Gamma"(link="log"), 
            data=finaldat, subset=Species_short=="AS")
summary(glm.AS)
Anova(glm.AS)
Tukey.testAS <- glht(glm.AS, linfct=mcp(Treatment = "Tukey")) 
summary(Tukey.testAS)

## Bidens pilosa
glm.BP<-glm(ExTotalB ~ Treatment + fNatB, family="Gamma"(link="log"), 
            data=finaldat, subset=Species_short=="BP")
summary(glm.BP)
Anova(glm.BP)
Tukey.testBP <- glht(glm.BP, linfct=mcp(Treatment="Tukey")) 
summary(Tukey.testBP)

## Conyza parva
glm.CP<-glm(ExTotalB ~ Treatment + fNatB, family="Gamma"(link="log"), 
            data=finaldat, subset=Species_short=="CP")
summary(glm.CP)
Anova(glm.CP)
Tukey.testCP <- glht(glm.CP, linfct=mcp(Treatment="Tukey")) 
summary(Tukey.testCP)

# Show model results in table ----------------------------------------
stargazer(glm.PA, glm.AS, glm.BP, glm.CP, align = TRUE, type = "html",
          out = "Output/Biomass_GLM_table.htm",
          dep.var.labels = "Exotic Plants Total Biomass",
          #covariate.labels = c("Salinity", "Native Saltmarsh Biomass", "Nutrients", "Intercept"),
          column.labels = c("P. aethiopicus", "A. sagittata", "B. bidens", "C. conyza")
          )


# Plotting of raw data ----------------------------------------------------
df.t <- ddply(finaldat, c( "Species","Treatment"), summarise,
              mean = mean(ExTotalB), median = median(ExTotalB), sd = sd(ExTotalB), 
              sem = sd(ExTotalB)/sqrt(length(ExTotalB)))
df.t[,3:6] <- round(df.t[,3:6], 2)

newRow1 <- data.frame(Species='Conyza parva', Treatment='HSHN', mean=0, median=0, sd=0, sem=0)
newRow2 <- data.frame(Species='Conyza parva', Treatment='HSNN', mean=0, median=0, sd=0, sem=0)
newRow3 <- data.frame(Species='Bidens pilosa', Treatment='HSNN', mean=0, median=0, sd=0, sem=0)
newRow4 <- data.frame(Species='Acetosa sagittata', Treatment='HSNN', mean=0, median=0, sd=0, sem=0)
df.t_plus <- rbind(df.t, newRow1, newRow2, newRow3, newRow4)

xaxistext <- element_text(face = "bold.italic", color = "black", size = 12)

t <- ggplot(df.t_plus, aes(Species, mean, fill = Treatment)) +
  geom_bar(stat= "identity", position=position_dodge(width=0.9), color = "black") +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0.2, position=position_dodge(width=0.9)) +
  scale_fill_grey(start = 1, end = 0.3) + scale_y_continuous(limit = c(0, 90)) +
  theme_classic(base_size = 16) +
  xlab("Exotic Species") + ylab("Total Biomass [mg]") + labs(fill="") +
  theme(axis.text.x = xaxistext, axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"), legend.text = element_text(size=14))

t
ggplot(df.t_plus, aes(Species, mean, fill = Treatment)) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0.2, position=position_dodge(width=0.5)) +
  geom_point(stat= "identity", position=position_dodge(width=0.5), shape = 21, size = 3, aes(fill=Treatment)) +
  scale_fill_grey(start = 1, end = 0.3) + scale_y_continuous(limit = c(0, 90)) +
  theme_classic(base_size = 16) +
  xlab("Exotic Species") + ylab("Total Biomass [mg]") + labs(fill="") +
  theme(axis.text.x = xaxistext, axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"), legend.text = element_text(size=14))



## add letters for significances

cld(Tukey.testPA)
cld(Tukey.testAS)
cld(Tukey.testBP)
cld(Tukey.testCP)

# explore glm data, incl. plotting, conf. intervals -----------------------
library(effects)

plot(glm.PA)
plot(allEffects(glm.PA),terms="Treatment", type="response", main="PA")
effects_glm.PA <- allEffects(glm.PA)
str(effects_glm.PA)
round(pR2(glm.PA),3) # pseudo R square for glm and round round to 3 digits
exp(effects_glm.PA$Treatment$lower) # lower confidence interval
exp(effects_glm.PA$Treatment$upper) # upper confidence interval

## originally, I wanted to do a Dunnett Test as post hoc test with HSNN as the control 
## group. However, this is not possible because there is not enough data per species.

# Set control group for Dunnett test, by default the first level would be chosen 
# as base line
control_group <- table(finaldat$Treatment) 
dunnett_grouping <- contrMat(control_group, base =2) #use 2 level as baseline
# GLM for PA
glm.PA<-glm(ExTotalB ~ Treatment + fNatB, family="Gamma"(link="log"), 
            data=finaldat, subset=Species_short=="PA")
summary(glm.PA)
Anova(glm.PA)
# Dunnett Test using the multcomp package and the second level in treatment as control group
Dunnett.testPA <- glht(glm.PA, linfct=mcp(Treatment=dunnett_grouping)) 
summary(Dunnett.testPA)