# Genotype ranking based on Clh A content results

# set working directory
setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/09.Manuscript/Stats_GCB_CBASSvsCLASSIC/Chlorophyll/")

#libraries to load
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)

# read working table
chlA.data <- read.delim("./Raw.data/Spis.CBASSvsCLASSIC.ChlA.data_17032019.txt", header = TRUE, sep = "\t")
colnames(chlA.data)

# keep only certain columns
chlA_data_long <- chlA.data %>% select(experiment, reef.site, temp.intensity, genotype, chlA_ug_cm2)

# change format from wide to long
chlA_data_wide <- spread(chlA_data_long, temp.intensity, chlA_ug_cm2)

# calculate the delta high
chlA_data_wide$deltahigh <- (chlA_data_wide$High - chlA_data_wide$Control)

# genotype correlation cbass vs classic - High ChlA content delta
high <- data.frame(genotype=chlA_data_wide$genotype[chlA_data_wide$experiment=="CBASS"],
                delta_high_x=chlA_data_wide$deltahigh[chlA_data_wide$experiment=="CBASS"],
                delta_high_y=chlA_data_wide$deltahigh[chlA_data_wide$experiment=="CLASSIC"],
                reef.site=chlA_data_wide$reef.site[chlA_data_wide$experiment=="CBASS"])

# correlation model - High ChlA content delta
cbass.vs.classixc_hightemp_corr <- cor.test(x = high$delta_high_x, y = high$delta_high_y, method = "pearson")
cbass.vs.classixc_hightemp_corr$p.value #0.8094154

cbass.vs.classixc_hightemp_lm <- lm(high$delta_high_x ~ high$delta_high_y)
summary(cbass.vs.classixc_hightemp_lm)
# Multiple R-squared:  0.005518
# p-value: 0.8094

# plot
Chlhigh <- ggplot(high,aes(x=delta_high_x,y=delta_high_y)) +
  theme_classic() +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_poly_eq(aes(label= paste(..eq.label..)), npcx = "right", npcy = 0.2, formula = y ~ x, parse=TRUE, size = 4)+
  stat_poly_eq(aes(label= paste(..rr.label..)), npcx = "right", npcy = 0.15, formula = y ~ x, parse=TRUE, size = 4)+
  stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),npcx = "right", npcy = 0.1, size = 4)+
  geom_point(aes(color=reef.site, shape=reef.site), size = 8, alpha=0.8)+
  theme(legend.position = 'bottom')+
  scale_x_continuous(limits = c(-2, 1.5))+
  scale_y_continuous(limits = c(-3, 0.1))+
  ggtitle("Chl. A content (HighTemp - ControlTemp)") + xlab("delta Chl. A (ug/cm2) CBASS") + ylab("delta Chl. A (ug/cm2) CLASSIC")+
  geom_text(aes(label=genotype), vjust=0, size=3)+
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=12, face="bold"),
        legend.text = element_text(colour="grey20", size=12, face="plain")) + scale_color_manual(values=c("#56B4E9", "#E69F00")) + scale_shape_manual(values=c(16, 16))

Chlhigh

ggsave("./genotype.ChlA/Plots/genotyperank_deltahighChlA_20200329.pdf", width = 6, height = 6)

