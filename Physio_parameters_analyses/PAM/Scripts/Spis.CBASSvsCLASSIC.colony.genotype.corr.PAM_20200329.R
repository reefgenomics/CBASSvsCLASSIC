# Genotype ranking based on PAM results

setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/09.Manuscript/Stats_GCB_CBASSvsCLASSIC/PAM/")

# libraries
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(svglite)

# data set
pam.data <- read.table("./Raw.data/CbassvsClassicMergedPAM_R_format_12082019.txt", header = TRUE)
colnames(pam.data)

# keep only certain columns
pam_data_long <- pam.data %>% select(experiment, reef.site, temp.intensity, replicate.group, genotype, yield)

# averaging PAm yield per genotype within each experiment (cbass or classic)
# average PAM yield per genptype within each experiment 
pam_average_geno <- aggregate(yield ~ reef.site + temp.intensity + experiment + genotype, data=pam_data_long, mean)

# change format from wide to long
pam_data_wide_genomean <- spread(pam_average_geno, temp.intensity, yield)

# calculate deltas
pam_data_wide_genomean$deltahigh <- (pam_data_wide_genomean$High - pam_data_wide_genomean$Control)
pam_data_wide_genomean$deltamedium <- (pam_data_wide_genomean$Medium - pam_data_wide_genomean$Control)


# plot genotype correlation by using a simple linear regression (square of Pearson's correlation coefficient is the same as the 𝑅2 in simple linear regression, as well as p-velue)

# cbass vs classic considering the average of replicate over genotype - High PAM delta # p-value added
high_genomean <- data.frame(genotype=pam_data_wide_genomean$genotype[pam_data_wide_genomean$experiment=="CBASS"], 
                            delta_high_x=pam_data_wide_genomean$deltahigh[pam_data_wide_genomean$experiment=="CBASS"], 
                            delta_high_y=pam_data_wide_genomean$deltahigh[pam_data_wide_genomean$experiment=="CLASSIC"], 
                            reef.site=pam_data_wide_genomean$reef.site[pam_data_wide_genomean$experiment=="CBASS"])

# correlation model - High PAM delta
cbass.vs.classixc_hightemp_corr <- cor.test(x = high_genomean$delta_high_x, y = high_genomean$delta_high_y, method = "pearson")
cbass.vs.classixc_hightemp_corr$p.value #0.004379476

cbass.vs.classixc_hightemp_lm <- lm(high_genomean$delta_high_x ~ high_genomean$delta_high_y)
cbass.vs.classixc_hightemp_lm <- lm(high_genomean$delta_high_y ~ high_genomean$delta_high_x)
summary(cbass.vs.classixc_hightemp_lm)
# Multiple R-squared:  0.537
# p-value: 0.004379

# plot
Phigh_genomean <- ggplot(high_genomean,aes(x=delta_high_x,y=delta_high_y)) +
  theme_classic() +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_poly_eq(aes(label= paste(..eq.label..)), npcx = "right", npcy = 0.2, formula = y ~ x, parse=TRUE, size = 4)+
  stat_poly_eq(aes(label= paste(..rr.label..)), npcx = "right", npcy = 0.15, formula = y ~ x, parse=TRUE, size = 4)+
  stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-value = ", signif(..p.value.., digits = 1), sep = "")),npcx = "right", npcy = 0.1, size = 4)+
  geom_point(aes(color=reef.site, shape=reef.site), size = 8, alpha=0.8)+
  theme(legend.position = 'bottom')+
  scale_x_continuous(limits = c(-0.2, 0.15))+
  scale_y_continuous(limits = c(-0.65, 0.1))+
  ggtitle("Fv/Fm (HighTemp - ControlTemp)") + xlab("delta Fv/Fm CBASS") + ylab("delta Fv/Fm CLASSIC")+
  geom_text(aes(label=genotype), vjust=0, size=3)+
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=20, face="bold"),
        legend.text = element_text(colour="grey20", size=18, face="plain")) + scale_color_manual(values=c("#56B4E9", "#E69F00")) + scale_shape_manual(values=c(16, 16))
Phigh_genomean

ggsave("./genotype.PAM/Plots/genotyperank_averagedgenotype_deltahighPAM_20200329.pdf", width = 6, height = 6)


# genotype correlation cbass vs classic considering considering the average of replicate over genotype - Medium PAM delta
medium_genomean <- data.frame(genotype=pam_data_wide_genomean$genotype[pam_data_wide_genomean$experiment=="CBASS"],
                  delta_medium_x=pam_data_wide_genomean$deltamedium[pam_data_wide_genomean$experiment=="CBASS"],
                  delta_medium_y=pam_data_wide_genomean$deltamedium[pam_data_wide_genomean$experiment=="CLASSIC"],
                  reef.site=pam_data_wide_genomean$reef.site[pam_data_wide_genomean$experiment=="CBASS"])

# correlation model - Medium PAM delta
cbass.vs.classixc_midtemp_corr <- cor.test(x = medium_genomean$delta_medium_x, y = medium_genomean$delta_medium_y, method = "pearson")
cbass.vs.classixc_midtemp_corr$p.value #0.9677844

cbass.vs.classixc_midtemp_lm <- lm(medium_genomean$delta_medium_x ~ medium_genomean$delta_medium_y)
summary(cbass.vs.classixc_midtemp_lm)
# Multiple R-squared:  0.0001417
# p-value: 0.9678

# plot
Pmedium_genomean=ggplot(medium_genomean,aes(x=delta_medium_x,y=delta_medium_y)) +
  theme_classic() +
  stat_smooth(method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_poly_eq(aes(label= paste(..eq.label..)), npcx = "right", npcy = 0.2, formula = y ~ x, parse=TRUE, size = 4)+
  stat_poly_eq(aes(label= paste(..rr.label..)), npcx = "right", npcy = 0.15, formula = y ~ x, parse=TRUE, size = 4)+
  stat_fit_glance(method = 'lm', method.args = list(formula=y ~ x), aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),npcx = "right", npcy = 0.1, size = 4)+
  geom_point(aes(color=reef.site, shape=reef.site), size = 8, alpha=0.8)+
  theme(legend.position = 'bottom')+
  scale_x_continuous(limits = c(-0.2, 0.15))+
  scale_y_continuous(limits = c(-0.65, 0.1))+
  ggtitle("Fv/Fm (MediumTemp - ControlTemp)") + xlab("delta Fv/Fm CBASS") + ylab("delta Fv/Fm CLASSIC")+
  geom_text(aes(label=genotype), vjust=0, size=3)+
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=20, face="bold"),
        legend.text = element_text(colour="grey20", size=18, face="plain"))+ scale_color_manual(values=c("#56B4E9", "#E69F00")) + scale_shape_manual(values=c(16, 16))
Pmedium_genomean

ggsave("./genotype.PAM/Plots/genotyperank_averagedgenotype_deltamediumPAM_20200329.pdf", width = 6, height = 6)


ggarrange(Pmedium_genomean, Phigh_genomean,
          ncol = 2, nrow = 1,
          common.legend = TRUE, legend = "bottom")

ggsave("./genotype.PAM/Plots/genotyperank_average_deltamedium_deltahighPAM_20200329.pdf", width = 12, height = 6)
