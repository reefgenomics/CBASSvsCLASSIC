# Stylophora pistillata - Tahla Reef
# Physiological measurements - Chlorophyll A

# Set working directory
setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/10.stats/Chlorophyll")

#libraries to load
library(ggplot2)
library(easyGgplot2)
library(dplyr)
library(car)
library(emmeans)

# Read working table
spis.tahla.physio.data <- read.delim("./Raw.data/Spis.CBASSvsCLASSIC.ChlA.data_17032019.txt", header = TRUE, sep = "\t")

# Set experiment and temperature intensity as factors
spis.tahla.physio.data$experiment_f <- factor(spis.tahla.physio.data$experiment, levels = c("CBASS", "CLASSIC"))
spis.tahla.physio.data$temp.intensity <-factor(spis.tahla.physio.data$temp.intensity, levels = c("Control", "Medium", "High", "Extreme"))

# CBASS data set
CBASS <- spis.tahla.physio.data %>%
  dplyr::filter(experiment == "CBASS")
# temperature intensity and reef site as factors
CBASS$temp.intensity <-factor(CBASS$temp.intensity, levels = c("Control", "Medium", "High", "Extreme"), ordered = TRUE)
CBASS$reef.site <-factor(CBASS$reef.site, levels = c("Exposed", "Protected"), ordered = TRUE)

# CLASSIC data set 
CLASSIC <- spis.tahla.physio.data %>%
  dplyr::filter(experiment == "CLASSIC")
# Temperature intensity and reef site as factors
CLASSIC$temp.intensity <-factor(CLASSIC$temp.intensity, levels = c("Control", "Medium", "High"), ordered = TRUE)
CLASSIC$reef.site <-factor(CLASSIC$reef.site, levels = c("Exposed", "Protected"), ordered = TRUE)

############CHLOROPHYLL CONTENT PER CM^2############
# This sample S4-E-ST-39B has a negative Chl A value, therefore it will be removed from the ChlA analysis (likely a technical mistake)
CBASS$chlA_ug_cm2[CBASS$chlA_ug_cm2<= 0] <- NA

########
# Statistical tests for CBASS chlA data
# Analysis of variance (2-way) ANOVA CBASS
chlA.CBASS.lm <- lm(chlA_ug_cm2 ~ reef.site * temp.intensity, data = CBASS)
Anova(chlA.CBASS.lm, type = "III")

# Model fitting and assumptions diagnostic 
plot(chlA_ug_cm2 ~ interaction(reef.site,temp.intensity), data = CBASS) # Box-plot homogeinity of variance
leveneTest(chlA_ug_cm2 ~ reef.site * temp.intensity, data=CBASS) # formal statistical test for homogeinity of variance
plot(chlA.CBASS.lm, 1) # Residual vs Fitted values
qqnorm(resid(chlA.CBASS.lm)); qqline(resid(chlA.CBASS.lm)) # qq plot to check for normal distribution of residuals
hist(resid(chlA.CBASS.lm)) # histogram of residuals to check for normal distribution of residuals
shapiro.test(chlA.CBASS.lm$residuals) # formal statistical test (not recommended due the small sample size)

###Check if by transforming the response variable, the residuals fit better to a normal distribution
# Square root transformation
chlA.CBASS.lm3 <- lm(sqrt(chlA_ug_cm2) ~ reef.site * temp.intensity, data = CBASS)
Anova(chlA.CBASS.lm, type = "III")

# Comparing Model fitting and assumptions diagnostic with and without sqrt transfromation
par(mfrow=c(1,2))
plot(chlA_ug_cm2 ~ interaction(reef.site,temp.intensity), data = CBASS) # Box-plot homogeinity of variance
plot(sqrt(chlA_ug_cm2) ~ interaction(reef.site,temp.intensity), data = CBASS)

par(mfrow=c(1,2))
plot(chlA.CBASS.lm, 1) # Residual vs Fitted values
plot(chlA.CBASS.lm3, 1)

par(mfrow=c(1,2))
qqnorm(resid(chlA.CBASS.lm)); qqline(resid(chlA.CBASS.lm)) # qq plot to check for normal distribution of residuals
qqnorm(resid(chlA.CBASS.lm3)); qqline(resid(chlA.CBASS.lm3))

par(mfrow=c(1,2))
hist(resid(chlA.CBASS.lm)) # histogram of residuals to check for normal distribution of residuals
hist(resid(chlA.CBASS.lm3))
dev.off()

# The square root transformation seems to help to the normal distribution of residuals and slightly to the homoscedasticity

# Post-hoc test - Meaningful comparisons
# Look for differences of reef sites within each temperature condition
chlA.CBASS.emms.reef <- emmeans(chlA.CBASS.lm3, pairwise ~ reef.site|temp.intensity, type="response", weights = "proportional", adjust="none")
summary(chlA.CBASS.emms.reef$emmeans)

# P.value adjustment of Bonferroni for the 9 aforemention comparisons. 
rbind(chlA.CBASS.emms.reef$contrasts, adjust="bonferroni")

#######
# Statistical tests for CLASSIC chlA data
# Analysis of variance (2-way) ANOVA CLASSIC
chlA.CLASSIC.lm <- lm(chlA_ug_cm2 ~ reef.site * temp.intensity, data = CLASSIC)
Anova(chlA.CLASSIC.lm, type = "III")

# Model fitting and assumptions diagnostic 
plot(chlA_ug_cm2 ~ interaction(reef.site,temp.intensity), data = CLASSIC) # Box-plot homogeinity of variance
leveneTest(chlA_ug_cm2 ~ reef.site * temp.intensity, data=CLASSIC) # formal statistical test for homogeinity of variance
plot(chlA.CLASSIC.lm, 1) # Residual vs Fitted values
qqnorm(resid(chlA.CLASSIC.lm)); qqline(resid(chlA.CLASSIC.lm)) # qq plot to check for normal distribution of residuals
hist(resid(chlA.CLASSIC.lm)) # histogram of residuals to check for normal distribution of residuals
shapiro.test(chlA.CLASSIC.lm$residuals) # formal statistical test (not recommended due the small sample size)

###Check if by transforming the response variable, the variance between treatments gets more homoscedastic
# Square root transformation
chlA.CLASSIC.lm3 <- aov(sqrt(chlA_ug_cm2) ~ reef.site * temp.intensity, data = CLASSIC)
Anova(chlA.CLASSIC.lm3, type = "III")

# Comparing Model fitting and assumptions diagnostic with and without sqrt transfromation
par(mfrow=c(1,2))
plot(chlA_ug_cm2 ~ interaction(reef.site,temp.intensity), data = CLASSIC) # Box-plot homogeinity of variance
plot(sqrt(chlA_ug_cm2) ~ interaction(reef.site,temp.intensity), data = CLASSIC)

par(mfrow=c(1,2))
plot(chlA.CLASSIC.lm, 1) # Residual vs Fitted values
plot(chlA.CLASSIC.lm3, 1)

par(mfrow=c(1,2))
qqnorm(resid(chlA.CLASSIC.lm)); qqline(resid(chlA.CLASSIC.lm)) # qq plot to check for normal distribution of residuals
qqnorm(resid(chlA.CLASSIC.lm3)); qqline(resid(chlA.CLASSIC.lm3))

par(mfrow=c(1,2))
hist(resid(chlA.CLASSIC.lm)) # histogram of residuals to check for normal distribution of residuals
hist(resid(chlA.CLASSIC.lm3))
dev.off()

# Post-hoc test - Meaningful comparisons
# Look for differences of reef sites within each temperature condition
chlA.CLASSIC.emms.reef <- emmeans(chlA.CLASSIC.lm3, pairwise ~ reef.site|temp.intensity, type="response", weights = "proportional", adjust="none")
summary(chlA.CLASSIC.emms.reef$emmeans)

# P.value adjustment of Bonferroni for the 9 aforemention comparisons. 
rbind(chlA.CLASSIC.emms.reef$contrasts, adjust="bonferroni")

# Box-Plot for Chl A content
chl <- ggplot(data=spis.tahla.physio.data, 
              aes(x=temp.intensity, y=chlA_ug_cm2, fill=reef.site)) +scale_fill_manual(values = c ("#56B4E9", "#E69F00"), name = "Reef site") +
  stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.7)+
  geom_boxplot(width=0.7, lwd=0.7, fatten=1) +
  expand_limits(y = 0)+
  facet_grid(~experiment_f, space = "free", scales = "free")+ 
  theme_bw()+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, colour = "black", size=2) +
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = 0, vjust = .5, face = "plain"),
        panel.spacing = unit(3, "lines"))

chl + xlab(label = "Temperature treatment") + ylab(label = "Chlorophyll A (ug/cm2)")+
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=20, face="bold"),
        legend.text = element_text(colour="grey20", size=18, face="plain"),
        legend.position="bottom")

ggsave("./Plots/CHLOROPH_CLASSIC.A_CBASS.B_22082019.pdf", width = 10, height = 6)
