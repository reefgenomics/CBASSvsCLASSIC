# Stylophora pistillata - Tahla Reef
# Physiological measurements - Symbiont density

# Set working directory
setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/10.stats/Cell.counts")

#libraries to load
library(ggplot2)
library(easyGgplot2)
library(dplyr)
library(car)
library(emmeans)

# Read working table
spis.tahla.physio.data <- read.delim("./Raw.data/Spis.CBASSvsCLASSIC.cellcounts.data_17032019.txt", header = TRUE, sep = "\t")

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

############SYMBIONT DENSITY############
# Statistical tests for CBASS symbiont density data
# Analysis of variance (2-way) ANOVA CBASS
symb.CBASS.lm <- lm(sym.cells_cm2 ~ reef.site * temp.intensity, data = CBASS)
Anova(symb.CBASS.lm, type = "III")

# Model fitting and assumptions diagnostic 
plot(sym.cells_cm2 ~ interaction(reef.site,temp.intensity), data = CBASS) # Box-plot homogeinity of variance
leveneTest(sym.cells_cm2 ~ reef.site * temp.intensity, data=CBASS) # formal statistical test for homogeinity of variance
plot(symb.CBASS.lm, 1) # Residual vs Fitted values
qqnorm(resid(symb.CBASS.lm)); qqline(resid(symb.CBASS.lm)) # qq plot to check for normal distribution of residuals
hist(resid(symb.CBASS.lm)) # histogram of residuals to check for normal distribution of residuals
shapiro.test(symb.CBASS.lm$residuals) # formal statistical test (not recommended due the small sample size)

# For each graphical plot of residuals, there is an associated formal statistical test which involves hypothesis testing (Seber, 1980).
# The main disadvantage of using a formal test is that sample size can largely affect the decision of whether the model fits the data or not (Cameron and Trivedi, 1998)
# According to many authors (e.g., Atkinson, 1987; Belsley, Kuh, & Welsch, 2005; Kozak, 2009; Moser & Stevens, 1992; Quinn & Keough, 2002; Rasch, Kubinger, & Moder, 2011; Schucany & Ng, 2006), significance tests should not be used for checking assumptions.
# Diagnostic residual plots are a better choice. On such diagnostic plots, we can identify untypical observations, heterogeneous within-treatment variances, and lack of normality of the residuals (and thus, of the dependent variable within treatments)

# A samples size of 7 is quite small, therefore it is advisable to based our assumption desition making on the residual plots instead of the stats

###Check if by transforming the response variable, the variance between treatments gets more homoscedastic
# Square root transformation
symb.CBASS.lm3 <- lm(sqrt(sym.cells_cm2) ~ reef.site * temp.intensity, data = CBASS)
Anova(symb.CBASS.lm3, type = "III")

# Comparing Model fitting and assumptions diagnostic with and without sqrt transfromation
par(mfrow=c(1,2))
plot(sym.cells_cm2 ~ interaction(reef.site,temp.intensity), data = CBASS) # Box-plot homogeinity of variance
plot(sqrt(sym.cells_cm2) ~ interaction(reef.site,temp.intensity), data = CBASS)

par(mfrow=c(1,2))
plot(symb.CBASS.lm, 1) # Residual vs Fitted values
plot(symb.CBASS.lm3, 1)

par(mfrow=c(1,2))
qqnorm(resid(symb.CBASS.lm)); qqline(resid(symb.CBASS.lm)) # qq plot to check for normal distribution of residuals
qqnorm(resid(symb.CBASS.lm3)); qqline(resid(symb.CBASS.lm3))

par(mfrow=c(1,2))
hist(resid(symb.CBASS.lm)) # histogram of residuals to check for normal distribution of residuals
hist(resid(symb.CBASS.lm3))
dev.off()

# The square root transformation seems to help to to the normal distribution of residuals and  homoscedasticity
# Two-way ANOVA test of the square root transformed data shows not significant differences between symbionts density means over temperature intensity, reef site, nor their interaction

# Post-hoc test - Meaningful comparisons
# Look for differences of reef sites within each temperature condition
#symb.cbass.emms.reef <- emmeans(symb.CBASS.lm3, pairwise ~ reef.site|temp.intensity, type="response")
#summary(symb.cbass.emms.reef$emmeans)

########
# Statistical tests for CLASSIC symbiont density data
# Analysis of variance (2-way) ANOVA CLASSIC
symb.CLASSIC.lm <- lm(sym.cells_cm2 ~ reef.site * temp.intensity, data = CLASSIC)
Anova(symb.CLASSIC.lm, type = "III")

# Model fitting and assumptions diagnostic 
plot(sym.cells_cm2 ~ interaction(reef.site,temp.intensity), data = CLASSIC) # Box-plot homogeinity of variance
leveneTest(sym.cells_cm2 ~ reef.site * temp.intensity, data=CLASSIC) # formal statistical test for homogeinity of variance
plot(symb.CLASSIC.lm, 1) # Residual vs Fitted values
qqnorm(resid(symb.CLASSIC.lm)); qqline(resid(symb.CLASSIC.lm)) # qq plot to check for normal distribution of residuals
hist(resid(symb.CLASSIC.lm)) # histogram of residuals to check for normal distribution of residuals
shapiro.test(symb.CLASSIC.lm$residuals) # formal statistical test (not recommended due the small sample size)

###Check if by transforming the response variable fits better to the normal distribution and the variance between treatments gets more homoscedastic
# Square root transformation
symb.CLASSIC.lm3 <- lm(sqrt(sym.cells_cm2) ~ reef.site * temp.intensity, data = CLASSIC)
Anova(symb.CLASSIC.lm3, type = "III")

# Comparing Model fitting and assumptions diagnostic with and without sqrt transfromation
par(mfrow=c(1,2))
plot(sym.cells_cm2 ~ interaction(reef.site,temp.intensity), data = CLASSIC) # Box-plot homogeinity of variance
plot(sqrt(sym.cells_cm2) ~ interaction(reef.site,temp.intensity), data = CLASSIC)

par(mfrow=c(1,2))
plot(symb.CLASSIC.lm, 1) # Residual vs Fitted values
plot(symb.CLASSIC.lm3, 1)

par(mfrow=c(1,2))
qqnorm(resid(symb.CLASSIC.lm)); qqline(resid(symb.CLASSIC.lm)) # qq plot to check for normal distribution of residuals
qqnorm(resid(symb.CLASSIC.lm3)); qqline(resid(symb.CLASSIC.lm3))

par(mfrow=c(1,2))
hist(resid(symb.CLASSIC.lm)) # histogram of residuals to check for normal distribution of residuals
hist(resid(symb.CLASSIC.lm3))
dev.off()

# The square root transformation seems to help to the normal distribution of residuals and slightly to the homoscedasticity
# There are significant differences (at the 0.05 level of significance) between the means of the symbiont density over temperature intensity and interaction between temperature intesity and reef

# Post-hoc test - Meaningful comparisons
# Look for differences of reef sites within each temperature condition
symb.classic.emms.reef <- emmeans(symb.CLASSIC.lm3, pairwise ~ reef.site|temp.intensity, type="response", weights = "proportional", adjust="none")
summary(symb.classic.emms.reef$emmeans)

# P.value adjustment of Bonferroni for the multiple comparisons. 
rbind(symb.classic.emms.reef$contrasts, adjust="bonferroni")


# Box-Plot for symbiont density
symb <- ggplot(data=spis.tahla.physio.data, 
               aes(x=temp.intensity, y=sym.cells_cm2, fill=reef.site)) +
  scale_fill_manual(values = c ("#56B4E9", "#E69F00"), name = "Reef site") +
  stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.7)+
  geom_boxplot(width=0.7, lwd=0.7, fatten=1) +
  expand_limits(y = c(0, 1200000))+
  facet_grid(~experiment_f, space = "free", scales = "free")+ #scales = free allows the removal of factor without data
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
        strip.text.x = element_text(color = "grey20", size = 9, angle = 0, hjust = 0, vjust = .5, face = "plain"),
        panel.spacing = unit(3, "lines"))

symb + xlab(label = "Temperature treatment") + ylab(label = "Symbiont density (cells/cm2)")+
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=20, face="bold"),
        legend.text = element_text(colour="grey20", size=18, face="plain"),
        legend.position="bottom")

ggsave("./Plots/SYMBdensit_CLASSIC.A_CBASS.B_22082019.pdf", width = 10, height = 6)
