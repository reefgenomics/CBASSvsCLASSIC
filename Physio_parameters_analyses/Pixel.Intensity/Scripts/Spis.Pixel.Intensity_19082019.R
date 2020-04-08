# Stylophora pistillata - Tahla Reef
# Pixel intensity analysis

setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/10.stats/Pixel.Intensity/")

# Libraries to load
library(ggplot2)
library(ggpubr)
library(easyGgplot2)
library(dplyr)
library(car)
library(emmeans)

# Read working table
spis.pixel <- read.table("./Raw.data/Spis.CLASSIC.A.CBASS.B_pixel.intensity_17032019.txt", header = TRUE, sep = "\t")

# Set experiment and temp.intensity as factors
spis.pixel$experiment <- factor(spis.pixel$experiment, levels = c("CBASS", "CLASSIC"), ordered = TRUE)
spis.pixel$temp.intensity <- factor(spis.pixel$temp.intensity, levels = c("Control", "Medium", "High", "Extreme"), ordered = TRUE)

# CBASS data set
CBASS <- spis.pixel %>%
  dplyr::filter(experiment == "CBASS")
# temperature intensity and reef site as factors
CBASS$temp.intensity <-factor(CBASS$temp.intensity, levels = c("Control", "Medium", "High", "Extreme"), ordered = TRUE)
CBASS$reef.site <-factor(CBASS$reef.site, levels = c("Exposed", "Protected"), ordered = TRUE)

# CLASSIC data set
CLASSIC <- spis.pixel %>%
  dplyr::filter(experiment == "CLASSIC")
# temperature intensity and reef site as factors
CLASSIC$temp.intensity <-factor(CLASSIC$temp.intensity, levels = c("Control", "Medium", "High"), ordered = TRUE)
CLASSIC$reef.site <-factor(CLASSIC$reef.site, levels = c("Exposed", "Protected"), ordered = TRUE)


############RED PIXELS############
# Statistical tests for CBASS RED pixel intensity data
# Analysis of variance (2-way) ANOVA CBASS
red.CBASS.lm <- lm(red ~ reef.site * temp.intensity, data = CBASS)
Anova(red.CBASS.lm, type = "III")

# Model fitting and assumptions diagnostic 
plot(red ~ interaction(reef.site,temp.intensity), data = CBASS) # Box-plot homogeinity of variance
leveneTest(red ~ reef.site * temp.intensity, data=CBASS) # formal statistical test for homogeinity of variance
plot(red.CBASS.lm, 1) # Residual vs Fitted values
qqnorm(resid(red.CBASS.lm)); qqline(resid(red.CBASS.lm)) # qq plot to check for normal distribution of residuals
hist(resid(pam.CBASS.lm)) # histogram of residuals to check for normal distribution of residuals
shapiro.test(red.CBASS.lm$residuals) # formal statistical test (not recommended due the small sample size)

# post hoc test - Meaningful comparisons
# Look for differences of reef sites within each temperature condition
red.CBASS.emms.reef <- emmeans(red.CBASS.lm, pairwise ~ reef.site|temp.intensity, type="response", weights = "proportional", adjust="none")
summary(red.CBASS.emms.reef$emmeans)
plot(red.CBASS.emms.reef)

# P.value adjustment of Bonferroni for the 9 aforemention comparisons. 
rbind(red.CBASS.emms.reef$contrasts, adjust="bonferroni")

# Stats CLASSIC RED-PIXEL
# Analysis of variance (2-way) ANOVA CLASSIC
red.CLASSIC.lm <- lm(red ~ reef.site * temp.intensity, data = CLASSIC)
Anova(red.CLASSIC.lm, type = "III")

# Model fitting and assumptions diagnostic 
plot(red ~ interaction(reef.site,temp.intensity), data = CLASSIC) # Box-plot homogeinity of variance
leveneTest(red ~ reef.site * temp.intensity, data=CLASSIC) # formal statistical test for homogeinity of variance
plot(red.CLASSIC.lm, 1) # Residual vs Fitted values
qqnorm(resid(red.CLASSIC.lm)); qqline(resid(red.CLASSIC.lm)) # qq plot to check for normal distribution of residuals
hist(resid(pam.CLASSIC.lm)) # histogram of residuals to check for normal distribution of residuals
shapiro.test(red.CLASSIC.lm$residuals) # formal statistical test (not recommended due the small sample size)

# post hoc test - Meaningful comparisons
# Look for differences of reef sites within each temperature condition
red.CLASSIC.emms.reef <- emmeans(red.CLASSIC.lm, pairwise ~ reef.site|temp.intensity, type="response", weights = "proportional", adjust="none")
summary(red.CLASSIC.emms.reef$emmeans)
plot(red.CLASSIC.emms.reef)

# P.value adjustment of Bonferroni for the 9 aforemention comparisons. 
rbind(red.CLASSIC.emms.reef$contrasts, adjust="bonferroni")

# box plot CBASS - CLASSIC analysis pixel red
spis.pix.red <- ggplot(data=spis.pixel, 
               aes(x=temp.intensity, y=red, fill=reef.site)) +
  scale_fill_manual(values = c ("#56B4E9", "#E69F00"), name = "Reef site") +
  stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.7)+
  geom_boxplot(width=0.7, lwd=0.7, fatten=1) +
  expand_limits(y = c(120, 255))+
  facet_grid(~experiment, space = "free", scales = "free")+ #scales = free allows the removal of factor without data
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
        strip.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = 0, vjust = .5, face = "plain"),
        panel.spacing = unit(3, "lines"))

spis.pix.red + xlab(label = "Temperature treatment") + ylab(label = "Intensity within the red channel\n(0 = darkest, 255 = brightest)")+
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=20, face="bold"),
        legend.text = element_text(colour="grey20", size=18, face="plain"),
        legend.position = "bottom")

ggsave("./Plots/Spis.PixelRed_CLASSIC.A_CBASS.B_22082019.pdf", width = 10, height = 6)
