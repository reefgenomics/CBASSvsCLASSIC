# Genotype ranking based on Red channel pixel intensity (Clh A content proxy) results

# set working directory
setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/09.Manuscript/Stats_GCB_CBASSvsCLASSIC/Pixel.Intensity/")

#libraries to load
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(svglite)

# read working table
Redpixel.data <- read.delim("./Raw.data/Spis.CLASSIC.A.CBASS.B_pixel.intensity_17032019.txt", header = TRUE, sep = "\t")
colnames(Redpixel.data)

# keep only certain columns
Redpixel_data_long <- Redpixel.data %>% select(experiment, reef.site, temp.intensity, genotype, red)

# change format from wide to long
Redpixel_data_wide <- spread(Redpixel_data_long, temp.intensity, red)

# calculate the delta high
Redpixel_data_wide$deltahigh <- (Redpixel_data_wide$High - Redpixel_data_wide$Control)

# genotype correlation cbass vs classic - High Red pixel delta
high <- data.frame(genotype=Redpixel_data_wide$genotype[Redpixel_data_wide$experiment=="CBASS"],
                   delta_high_x=Redpixel_data_wide$deltahigh[Redpixel_data_wide$experiment=="CBASS"],
                   delta_high_y=Redpixel_data_wide$deltahigh[Redpixel_data_wide$experiment=="CLASSIC"],
                   reef.site=Redpixel_data_wide$reef.site[Redpixel_data_wide$experiment=="CBASS"])

# correlation model - High Red pixel intensity delta
cbass.vs.classixc_hightemp_corr <- cor.test(x = high$delta_high_x, y = high$delta_high_y, method = "pearson")
cbass.vs.classixc_hightemp_corr$p.value #0.3570165

cbass.vs.classixc_hightemp_lm <- lm(high$delta_high_x ~ high$delta_high_y)
summary(cbass.vs.classixc_hightemp_lm)
# Multiple R-squared:  0.07751
# p-value: 0.357

# plot
Redpixelhigh <- ggplot(high,aes(x=delta_high_x,y=delta_high_y)) +
  theme_classic() +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_poly_eq(aes(label= paste(..eq.label..)), npcx = "right", npcy = 0.2, formula = y ~ x, parse=TRUE, size = 4)+
  stat_poly_eq(aes(label= paste(..rr.label..)), npcx = "right", npcy = 0.15, formula = y ~ x, parse=TRUE, size = 4)+
  stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),npcx = "right", npcy = 0.1, size = 4, family="Helvetica", colour="grey20")+
  geom_point(aes(color=reef.site, shape=reef.site), size = 8, alpha=0.8)+
  theme(legend.position = 'bottom')+
  ggtitle("Red pixel intensity (HighTemp - ControlTemp)") + xlab("delta Intensity Red Channel CBASS") + ylab("delta Intensity Red Channel CLASSIC")+
  geom_text(aes(label=genotype), vjust=0, size=3)+
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.ticks.length=unit(.35, "cm"),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain", family="Helvetica"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain", family="Helvetica"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain", family="Helvetica"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain", family="Helvetica"),
        legend.title = element_text(colour="grey20", size=12, face="bold"),
        legend.text = element_text(colour="grey20", size=12, face="plain")) + scale_color_manual(values=c("#56B4E9", "#E69F00")) + scale_shape_manual(values=c(16, 16))

Redpixelhigh

ggsave("./genotype.Pixel/Plots/genotyperank_deltahighRedPixel_20200329.pdf", width = 6, height = 6)


