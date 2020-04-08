# Stylophora pistillata - Tahla Reef
# Physiological measurements - PAM

setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/09.Manuscript/Stats_GCB_CBASSvsCLASSIC/PAM/")

#libraries
library(car)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(sjPlot)
library(sjmisc)

# data set
pam.data <- read.table("./Raw.data/CbassvsClassicMergedPAM_R_format_12082019.txt", header = T)
colnames(pam.data)

pam.data$genotype <- as.factor(pam.data$genotype)
pam.data$reef.site <-  as.factor(pam.data$reef.site)
pam.data$replicate.group <- as.factor(pam.data$replicate.group)

########
# Statistical tests for CBASS PAM data
# Analysis of variance (2-way) ANOVA CBASS (mixed linear model)
cbass <- subset(pam.data, pam.data$experiment == "CBASS")
cbass$temperature <- as.factor(cbass$temperature)

# Mixed Effects model - random effect of replicate.group
# Since the model is completely balanced the SS type I, II or III will provide the same result
pam.cbass.model <- lmer(yield ~ reef.site * temperature + (1|tank), data = cbass)
anova(pam.cbass.model) # anova from lmerTest ANOVA SS Type III with ddf="Satterthwaite"
ranova(pam.cbass.model) # anova-like table for random effects

# Model fitting and assumptions diagnostic 
plot(yield ~ interaction(reef.site,temperature,replicate.group), data = cbass) # Box-plot homogeinity of variance
leveneTest(yield ~ reef.site * temperature * replicate.group, data=cbass) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
plot(pam.cbass.model) # Residual vs Fitted values
qqnorm(resid(pam.cbass.model)); qqline(resid(pam.cbass.model)) # qq plot to check for normal distribution of residuals
hist(resid(pam.cbass.model)) # histogram of residuals to check for normal distribution of residuals
shapiro.test(resid(pam.cbass.model)) # formal statistical test (not recommended due the small sample size)

# comparing between reef sites within each temperature treatment
pam.cbass.emms.reef <- emmeans(pam.cbass.model, pairwise ~ reef.site|temperature, weights = "proportional", adjust="none")
summary(pam.cbass.emms.reef$emmeans)

# P.value adjustment of the Bonferroni
rbind(pam.cbass.emms.reef$contrasts, adjust="bonferroni")

########
# Statistical tests for CLASSIC symbiont density data
# Analysis of variance (2-way) ANOVA CLASSIC (mixed linear model)
classic <- subset(pam.data, pam.data$experiment == "CLASSIC")
classic$temperature <- as.factor(classic$temperature)

# Mixed Effects model - random effect of replicate.group
pam.classic.model <- lmer(yield ~ reef.site * temperature + (1 | tank), data = classic)
anova(pam.classic.model) # anova from lmerTest ANOVA SS Type III with ddf="Satterthwaite"
ranova(pam.classic.model) # anova-like table for random effects

# Model fitting and assumptions diagnostic 
plot(yield ~ interaction(temperature,reef.site,replicate.group), data = classic) # Box-plot homogeinity of variance
leveneTest(yield ~ reef.site * temperature * replicate.group, data=classic) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
plot(pam.classic.model) # Residual vs Fitted values
qqnorm(resid(pam.classic.model)); qqline(resid(pam.classic.model)) # qq plot to check for normal distribution of residuals
hist(resid(pam.classic.model)) # histogram of residuals to check for normal distribution of residuals
shapiro.test(resid(pam.classic.model)) # formal statistical test (not recommended due the small sample size)

# comparing between reef sites within each temperature treatment
pam.classic.emms.reef <- emmeans(pam.classic.model, pairwise ~ reef.site|temperature, weights = "proportional", adjust="none")
summary(pam.classic.emms.reef$emmeans)

# P.value adjustment of Bonferroni
rbind(pam.classic.emms.reef$contrasts, adjust="bonferroni")

### BOX plots
pam.data$experiment <- as.factor(pam.data$experiment)
pam.data$temp.intensity <- as.factor(pam.data$temp.intensity)
print(levels(pam.data$temp.intensity))
pam.data$temp.intensity = factor(pam.data$temp.intensity,levels(pam.data$temp.intensity)[c(1,4,3:2)]) # reorder the level of a factor

pam <- ggplot(data=pam.data, 
              aes(x=temp.intensity, y=yield, label= temperature, fill=reef.site)) +
  scale_fill_manual(values = c ("#56B4E9", "#E69F00"), name = "Reef site") +
  stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.7)+
  geom_boxplot(width=0.7, lwd=0.7, fatten=1) +
  expand_limits(y = 0)+
  facet_grid(~experiment, space = "free", scales = "free")+ #this can also be used but rows and columns have to be specified -> facet_grid(. ~ experiment)
  theme_bw() +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, colour = "black", size=2) +
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = 0, vjust = 0.5, face = "plain"),
        panel.spacing = unit(3, "lines"))
  

pam + xlab(label = "Temperature treatment") + ylab(label = "Photosynthetic efficiency (Fv/Fm)")+
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=20, face="bold"),
        legend.text = element_text(colour="grey20", size=18, face="plain"),
        legend.position="bottom")

ggsave("./Plots/PAMdata_CLASSIC_CBASS_22082019.pdf", width = 10, height = 6)
