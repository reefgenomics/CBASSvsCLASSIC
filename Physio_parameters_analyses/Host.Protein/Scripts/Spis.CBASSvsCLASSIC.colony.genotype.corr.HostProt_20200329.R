# Genotype ranking based on Host Protein results

# set working directory
setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/09.Manuscript/Stats_GCB_CBASSvsCLASSIC/Host.Protein/")

# libraries to load
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)

# data set
Hostprot.data <- read.delim("./Raw.data/Spis.CBASSvsCLASSIC.HostProtein.data_17032019.txt", header = TRUE, sep = "\t")
colnames(Hostprot.data)

# keep only certain columns
Hostprot_data_long <- Hostprot.data %>% select(experiment, reef.site, temp.intensity, genotype, host.protein_mg_cm2)

# change format from wide to long
Hostprot_data_wide <- spread(Hostprot_data_long, temp.intensity, host.protein_mg_cm2)

# calculate the delta high
Hostprot_data_wide$deltahigh <- (Hostprot_data_wide$High - Hostprot_data_wide$Control)

# genotype correlation cbass vs classic - High host protein delta
high <- data.frame(genotype=Hostprot_data_wide$genotype[Hostprot_data_wide$experiment=="CBASS"],
                delta_high_x=Hostprot_data_wide$deltahigh[Hostprot_data_wide$experiment=="CBASS"],
                delta_high_y=Hostprot_data_wide$deltahigh[Hostprot_data_wide$experiment=="CLASSIC"],
                reef.site=Hostprot_data_wide$reef.site[Hostprot_data_wide$experiment=="CBASS"])

# correlation model - High host protein delta
cbass.vs.classixc_hightemp_corr <- cor.test(x = high$delta_high_x, y = high$delta_high_y, method = "pearson")
cbass.vs.classixc_hightemp_corr$p.value #0.5496175

cbass.vs.classixc_hightemp_lm <- lm(high$delta_high_x ~ high$delta_high_y)
summary(cbass.vs.classixc_hightemp_lm)

# plot
HPhigh <- ggplot(high,aes(x=delta_high_x,y=delta_high_y)) +
  theme_classic() +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_poly_eq(aes(label= paste(..eq.label..)), npcx = "right", npcy = 0.2, formula = y ~ x, parse=TRUE, size = 4)+
  stat_poly_eq(aes(label= paste(..rr.label..)), npcx = "right", npcy = 0.15, formula = y ~ x, parse=TRUE, size = 4)+
  stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), npcx = "right", npcy = 0.1, size = 4)+
  geom_point(aes(color=reef.site, shape=reef.site), size = 8, alpha=0.8)+
  theme(legend.position = 'bottom')+
  scale_x_continuous(limits = c(-0.24, 0.16))+
  scale_y_continuous(limits = c(-0.24, 0.2))+
  ggtitle("Host protein content (HighTemp - ControlTemp)") + xlab("delta Host protein (mg/cm2) CBASS") + ylab("delta Host protein (mg/cm2) CLASSIC")+
  geom_text(aes(label=genotype), vjust=0, size=3)+
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=12, face="bold"),
        legend.text = element_text(colour="grey20", size=12, face="plain")) + scale_color_manual(values=c("#56B4E9", "#E69F00")) + scale_shape_manual(values=c(16, 16))

HPhigh

ggsave("./genotype.HostProt/Plots/genotyperank_deltahighHostprot_20200329.pdf", width = 6, height = 6)



