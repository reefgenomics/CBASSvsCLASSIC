# Genotype ranking based on Symbiodiniacea density results

# set working directory
setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/09.Manuscript/Stats_GCB_CBASSvsCLASSIC/Cell.counts/")

#libraries to load
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)

# read working table
cellcount.data <- read.delim("./Raw.data/Spis.CBASSvsCLASSIC.cellcounts.data_17032019.txt", header = TRUE, sep = "\t")
colnames(cellcount.data)

# keep only certain columns
cellcount_data_long <- cellcount.data %>% select(experiment, reef.site, temp.intensity, genotype, sym.cells_cm2)

# change format from wide to long
cellcount_data_wide <- spread(cellcount_data_long, temp.intensity, sym.cells_cm2)

# calculate the delta high
cellcount_data_wide$deltahigh <- (cellcount_data_wide$High - cellcount_data_wide$Control)

# genotype correlation cbass vs classic - High symbiont density delta
high <- data.frame(genotype=cellcount_data_wide$genotype[cellcount_data_wide$experiment=="CBASS"],
                   delta_high_x=cellcount_data_wide$deltahigh[cellcount_data_wide$experiment=="CBASS"],
                   delta_high_y=cellcount_data_wide$deltahigh[cellcount_data_wide$experiment=="CLASSIC"],
                   reef.site=cellcount_data_wide$reef.site[cellcount_data_wide$experiment=="CBASS"])

# correlation model - High symbiont density delta
cbass.vs.classixc_hightemp_corr <- cor.test(x = high$delta_high_x, y = high$delta_high_y, method = "pearson")
cbass.vs.classixc_hightemp_corr$p.value #0.9276077

cbass.vs.classixc_hightemp_lm <- lm(high$delta_high_x ~ high$delta_high_y)
summary(cbass.vs.classixc_hightemp_lm)
# Multiple R-squared:  0.000785
# p-value: 0.9276

# plot
Cellcounthigh <- ggplot(high,aes(x=delta_high_x,y=delta_high_y)) +
  theme_classic() +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_poly_eq(aes(label= paste(..eq.label..)), npcx = "right", npcy = 0.2, formula = y ~ x, parse=TRUE, size = 4)+
  stat_poly_eq(aes(label= paste(..rr.label..)), npcx = "right", npcy = 0.15, formula = y ~ x, parse=TRUE, size = 4)+
  stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),npcx = "right", npcy = 0.1, size = 4)+
  geom_point(aes(color=reef.site, shape=reef.site), size = 8, alpha=0.8)+
  theme(legend.position = 'bottom')+
  ggtitle("Symbiodinaceae density (HighTemp - ControlTemp)") + xlab("delta Symb. density (cells/cm2) CBASS") + ylab("delta Symb. density (cells/cm2) CLASSIC")+
  geom_text(aes(label=genotype), vjust=0, size=3)+
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=12, face="bold"),
        legend.text = element_text(colour="grey20", size=12, face="plain")) + scale_color_manual(values=c("#56B4E9", "#E69F00")) + scale_shape_manual(values=c(16, 16))

Cellcounthigh

ggsave("./genotype.Cell.counts/Plots/genotyperank_deltahighCellcount_20200329.pdf", width = 6, height = 6)
