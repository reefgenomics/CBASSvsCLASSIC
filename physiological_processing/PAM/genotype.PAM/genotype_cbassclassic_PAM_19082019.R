# Genotype ranking based on PAM results

setwd("~/Documents/Barshis-project/03.SUMMER-CRUISE/04.Short-Long.term.heat.stress.experiment/10.stats/PAM")

#libraries
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# data set
pam.data <- read.table("./Raw.data/CbassvsClassicMergedPAM_R_format_12082019.txt", header = T)
colnames(pam.data)

# keep only certain columns
pam_data_long <- pam.data %>% select(experiment, reef.site, temp.intensity, replicate.group, genotype, yield)

# Genotype corelation with replicate A and B separadtely
# Change format from wide to long
pam_data_wide <- spread(pam_data_long, temp.intensity, yield)

pam_data_wide$genorep <- paste(pam_data_wide$genotype, pam_data_wide$replicate.group, sep = "_")

pam_data_wide$deltahigh <- (pam_data_wide$High - pam_data_wide$Control)
pam_data_wide$deltamedium <- (pam_data_wide$Medium - pam_data_wide$Control)

# Plot genotype correlation cbass vs classic considering both replicates - High PAM delta
high=data.frame(genotype=pam_data_wide$genotype[pam_data_wide$experiment=="CBASS"],
                delta_high_x=pam_data_wide$deltahigh[pam_data_wide$experiment=="CBASS"],
                delta_high_y=pam_data_wide$deltahigh[pam_data_wide$experiment=="CLASSIC"],
                reef.site=pam_data_wide$reef.site[pam_data_wide$experiment=="CBASS"],
                genorep=pam_data_wide$genorep[pam_data_wide$experiment=="CBASS"])

Phigh=ggplot(high,aes(x=delta_high_x,y=delta_high_y)) +
  theme_classic() +
  stat_smooth(geom = "smooth", position = "identity", method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")), formula = y ~ x) +
  geom_point(aes(color=reef.site, shape=reef.site), size = 8, alpha=0.8)+
  theme(legend.position = 'bottom')+
  scale_x_continuous(limits = c(-0.2, 0.15))+
  scale_y_continuous(limits = c(-0.65, 0.1))+
  ggtitle("Fv/Fm (HighTemp - ControlTemp)") + xlab("delta Fv/Fm CBASS") + ylab("delta Fv/Fm CLASSIC")+
  geom_text(aes(label=genorep), vjust=-1.5, size=3)+
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=12, face="bold"),
        legend.text = element_text(colour="grey20", size=12, face="plain")) + scale_color_manual(values=c("#56B4E9", "#E69F00")) + scale_shape_manual(values=c(16, 16))

Phigh

ggsave("./genotype.PAM/Plots/genotyperank_bothreplicates_deltahighPAM_20180826.pdf", width = 6, height = 6)

# Plot genotype correlation cbass vs classic considering both replicates - Medium PAM delta
medium=data.frame(genotype=pam_data_wide$genotype[pam_data_wide$experiment=="CBASS"],
                delta_medium_x=pam_data_wide$deltamedium[pam_data_wide$experiment=="CBASS"],
                delta_medium_y=pam_data_wide$deltamedium[pam_data_wide$experiment=="CLASSIC"],
                reef.site=pam_data_wide$reef.site[pam_data_wide$experiment=="CBASS"],
                genorep=pam_data_wide$genorep[pam_data_wide$experiment=="CBASS"])
Pmedium=ggplot(medium,aes(x=delta_medium_x,y=delta_medium_y)) +
  theme_classic() +
  stat_smooth(geom = "smooth",position = "identity", method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")), formula = y ~ x) +
  geom_point(aes(color=reef.site, shape=reef.site), size = 8, alpha=0.8)+
  theme(legend.position = 'bottom')+
  scale_x_continuous(limits = c(-0.2, 0.15))+
  scale_y_continuous(limits = c(-0.65, 0.1))+
  ggtitle("Fv/Fm (MediumTemp - ControlTemp)") + xlab("delta Fv/Fm CBASS") + ylab("delta Fv/Fm CLASSIC")+
  geom_text(aes(label=genorep), vjust=-1.5, size=3)+
  theme(line= element_line(size = 1),
        axis.line = element_line(colour = "grey20"),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=12, face="bold"),
        legend.text = element_text(colour="grey20", size=12, face="plain"))+ scale_color_manual(values=c("#56B4E9", "#E69F00")) + scale_shape_manual(values=c(16, 16))

Pmedium 

ggsave("./genotype.PAM/Plots/genotyperank_bothreplicates_deltamediumPAM_20180826.pdf", width = 6, height = 6)

#####################################################################
# Averaging PAm yield per genotype within each experiment (cbass or classic)
# Average PAM yield per genptype within each experiment 
pam_average_geno <- aggregate(yield ~ reef.site + temp.intensity + experiment + genotype, data=pam_data_long, mean)

# Change format from wide to long
pam_data_wide_genomean <- spread(pam_average_geno, temp.intensity, yield)

pam_data_wide_genomean$genorep <- paste(pam_data_wide_genomean$genotype, pam_data_wide_genomean$replicate.group, sep = "_")

pam_data_wide_genomean$deltahigh <- (pam_data_wide_genomean$High - pam_data_wide_genomean$Control)
pam_data_wide_genomean$deltamedium <- (pam_data_wide_genomean$Medium - pam_data_wide_genomean$Control)

# Plot genotype correlation cbass vs classic considering the average of replicate over genotype - High PAM delta
high_genomean=data.frame(genotype=pam_data_wide_genomean$genotype[pam_data_wide_genomean$experiment=="CBASS"],
                delta_high_x=pam_data_wide_genomean$deltahigh[pam_data_wide_genomean$experiment=="CBASS"],
                delta_high_y=pam_data_wide_genomean$deltahigh[pam_data_wide_genomean$experiment=="CLASSIC"],
                reef.site=pam_data_wide_genomean$reef.site[pam_data_wide_genomean$experiment=="CBASS"])

Phigh_genomean=ggplot(high_genomean,aes(x=delta_high_x,y=delta_high_y)) +
  theme_classic() +
  stat_smooth(geom = "smooth",position = "identity", method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")), formula = y ~ x) +
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

ggsave("./genotype.PAM/Plots/genotyperank_averagedgenotype_deltahighPAM_20180826().pdf", width = 6, height = 6)

# Plot genotype correlation cbass vs classic considering considering the average of replicate over genotype - Medium PAM delta
medium_genomean=data.frame(genotype=pam_data_wide_genomean$genotype[pam_data_wide_genomean$experiment=="CBASS"],
                  delta_medium_x=pam_data_wide_genomean$deltamedium[pam_data_wide_genomean$experiment=="CBASS"],
                  delta_medium_y=pam_data_wide_genomean$deltamedium[pam_data_wide_genomean$experiment=="CLASSIC"],
                  reef.site=pam_data_wide_genomean$reef.site[pam_data_wide_genomean$experiment=="CBASS"])

Pmedium_genomean=ggplot(medium_genomean,aes(x=delta_medium_x,y=delta_medium_y)) +
  theme_classic() +
  stat_smooth(geom = "smooth",position = "identity", method = "lm",formula = y ~ x, se = TRUE, level = 0.95, na.rm = TRUE, colour="grey40")+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")), formula = y ~ x) +
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
ggsave("./genotype.PAM/Plots/genotyperank_averagedgenotype_deltamediumPAM_20180826().pdf", width = 6, height = 6)


####

ggarrange(Pmedium, Phigh,
          ncol = 2, nrow = 1,
          common.legend = TRUE, legend = "bottom")

ggsave("./genotype.PAM/Plots/genotyperank_replicates_deltamedium_deltahighPAM_20180826.pdf", width = 12, height = 6)

ggarrange(Pmedium_genomean, Phigh_genomean,
          ncol = 2, nrow = 1,
          common.legend = TRUE, legend = "bottom")

ggsave("./genotype.PAM/Plots/genotyperank_average_deltamedium_deltahighPAM_20180826().pdf", width = 12, height = 6)
