# Seedling trait shifts of Sandberg bluegrass (POSE) grown in the greenhouse
# Two treatments: Cheatgrass(BRTE) competition and water availability

# Set the data pathway first

# Analysis:
# Did POSE traits shift by BRTE competition and water treatment?
# How did each trait respond to each treatment?
# Which trait is correlated with higher POSE performance?
# Is trait plasticity correlated with drought tolerance?

# Data
source("data_compiling/compile_demography.R")

# Packages
library(tidyverse) #data wrangling
library(ggplot2) #plot
library(ggpubr) #combine plots
library(vegan) #pca
library(corrplot) #correlation matrix
library(dplyr)
library(codyn)

# Prep data
root_traits <- full_join(rootbiomass, root) %>%
  group_by(PotID) %>%
  mutate(SRL = Length_cm/Root_Weight_g,
         Coarse = LenTotHistoClasses-D0L05,
         Fine = D0L05,
         PropF = D0L05/LenTotHistoClasses) %>% #calculate specific root length, coarse root length, and fine root length
  dplyr::select(PotID, Length_cm, SRL, SurfArea_cm2, AvgDiam_mm, Tips, Forks, Fine, Coarse, PropF) 
  
biomass_traits <- biomass %>%
  filter(Species == "POSE") %>%
  inner_join(., rootbiomass) %>%
  dplyr::select(-Species) %>%
  mutate(TotalBiomass = Dry_Biomass_Weight_g+Root_Weight_g,
         RMR = Root_Weight_g/TotalBiomass) %>% #calculate total biomass weight and root to mass ratio
  dplyr::select(PotID, TotalBiomass, RMR)

height_trait<- height %>% 
  group_by(PotID) %>% 
  summarise(Height = mean(Height_cm)) #calculate average seedling height

leaf_traits <- leaftraits %>%
  drop_na() %>%
  group_by(PotID) %>%
  summarise(SLA = FreshLeafArea_cm2/DryLeafWeight_g/NumberLeaves, 
         LDMC = DryLeafWeight_g/FreshLeafArea_cm2/NumberLeaves) #calculate specific leaf area and leaf dry matter content

germination_dates <- germination %>%
  dplyr::select(PotID, Days_emergence) %>%
  drop_na()

trait_master <- inner_join(potID, root_traits)%>%
  inner_join(.,biomass_traits) %>%
  inner_join(., height_trait) %>%
  inner_join(., leaf_traits) %>%
  inner_join(., germination_dates) %>%
  filter(Population != "Gund")
colnames(trait_master) <- c('PotID', 'Life_stage', 'Competition', 'Water', 'Population', 'Replicate',
                            'Length', 'SRL', 'SurfArea', 'AvgDiam', 'Tips', 'Forks', 'Fine', 'Coarse', 'PropF',
                            'TotalBiomass', 'RMR', 'Height', 'SLA', 'LDMC', 'Emergence')

trait_matrix_raw <- as.matrix(trait_master[,7:ncol(trait_master)]) #extract only the trait value columns

growth <- inner_join(potID, stemcount) %>% filter(Population != "Gund")

biomass_dat <- inner_join(potID, biomass) %>%
  filter(Population != "Gund") %>%
  pivot_wider(names_from = Species, values_from = Dry_Biomass_Weight_g)

# Function for standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
} 

#-----------------------------------------
# Take out any traits that are products of size and keep relativized metrics
# Standardize the data
pairs(~RMR + Tips + Length + Fine + Coarse + SRL + SurfArea + AvgDiam + Forks +PropF + TotalBiomass, trait_matrix_raw)
trait_matrix <- as.data.frame(decostand(trait_matrix_raw, "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)

#----------------------------------------------------------#
# Did POSE traits shift by BRTE competition and water treatment?
# Create PCA of traits
pca_trait = rda(trait_matrix, scale = FALSE) #run PCA on all traits
biplot(pca_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores <- as.data.frame(scores(pca_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab = as.data.frame(cbind(trait_master[,1:6],pca_trait_scores)) #add plot info back
pca_trait_scores_lab$Population <- ordered(as.factor(pca_trait_scores_lab$Population), levels = c("Butte Valley","Steens","EOARC",
                                                                                                  "Water Canyon",  "Reno"))
pca_trait_scores_lab$Treatment <- apply(pca_trait_scores_lab[ ,3:4 ] , 1 , paste , collapse = "_" )
envout<-as.data.frame(scores(pca_trait, choices=c(1,2), display=c("species")))
summary(pca_trait)
ggplot(pca_trait_scores_lab, aes(x = PC1, y = PC2))+
  geom_point(size = 4, aes(col = Treatment, shape = Population), alpha = 0.5)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15))+
  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
            fontface = "bold", label = row.names(envout), size = 5)+
  xlim(-2, 2.3)+
  xlab("PC1 (39.8%)")+
  ylab("PC2 (16.0%)")+
  scale_color_manual(values=c("#999999", "#6A0DAD", "#E69F00", "#56B4E9"))

# PCA within each population:
# Butte Valley
trait_Butte <- trait_master %>%
  filter(Population == "Butte Valley") %>%
  as.data.frame(decostand(., "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
trait_Butte_matrix <- as.matrix(trait_Butte[,7:ncol(trait_Butte)])
pca_Butte_trait = rda(trait_Butte_matrix , scale = TRUE) #run PCA on all traits
biplot(pca_Butte_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores_Butte <- as.data.frame(scores(pca_Butte_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab_Butte = as.data.frame(cbind(trait_Butte[,1:6],pca_trait_scores_Butte)) #add plot info back
pca_trait_scores_lab_Butte$Treatment <- apply(pca_trait_scores_lab_Butte[ ,3:4 ] , 1 , paste , collapse = "_" )
envout<-as.data.frame(scores(pca_Butte_trait , choices=c(1,2), display=c("species")))
summary(pca_Butte_trait)
PCA_Butte <- ggplot(pca_trait_scores_lab_Butte, aes(x = PC1, y = PC2))+
  geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15))+
  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
            fontface = "bold", label = row.names(envout), size = 5)+
  #xlim(-2, 2.3)+
  xlab("PC1 (32.7%)")+
  ylab("PC2 (21.9%)")

# Steens
trait_Steens <- trait_master %>%
  filter(Population == "Steens") %>%
  as.data.frame(decostand(., "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
trait_Steens_matrix <- as.matrix(trait_Steens[,7:ncol(trait_Steens)])
pca_Steens_trait = rda(trait_Steens_matrix , scale = TRUE) #run PCA on all traits
biplot(pca_Steens_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores_Steens <- as.data.frame(scores(pca_Steens_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab_Steens = as.data.frame(cbind(trait_Steens[,1:6],pca_trait_scores_Steens)) #add plot info back
pca_trait_scores_lab_Steens$Treatment <- apply(pca_trait_scores_lab_Steens[ ,3:4 ] , 1 , paste , collapse = "_" )
envout<-as.data.frame(scores(pca_Steens_trait , choices=c(1,2), display=c("species")))
summary(pca_Steens_trait)
PCA_Steens <- ggplot(pca_trait_scores_lab_Steens, aes(x = PC1, y = PC2))+
  geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15))+
  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
            fontface = "bold", label = row.names(envout), size = 5)+
  #xlim(-2, 2.3)+
  xlab("PC1 (42.1%)")+
  ylab("PC2 (18.4%)")

# EOARC
trait_EOARC <- trait_master %>%
  filter(Population == "EOARC") %>%
  as.data.frame(decostand(., "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
trait_EOARC_matrix <- as.matrix(trait_EOARC[,7:ncol(trait_EOARC)])
pca_EOARC_trait = rda(trait_EOARC_matrix , scale = TRUE) #run PCA on all traits
biplot(pca_EOARC_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores_EOARC <- as.data.frame(scores(pca_EOARC_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab_EOARC = as.data.frame(cbind(trait_EOARC[,1:6],pca_trait_scores_EOARC)) #add plot info back
pca_trait_scores_lab_EOARC$Treatment <- apply(pca_trait_scores_lab_EOARC[ ,3:4 ] , 1 , paste , collapse = "_" )
envout<-as.data.frame(scores(pca_EOARC_trait , choices=c(1,2), display=c("species")))
summary(pca_EOARC_trait)
PCA_EOARC <- ggplot(pca_trait_scores_lab_EOARC, aes(x = PC1, y = PC2))+
  geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15))+
  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
            fontface = "bold", label = row.names(envout), size = 5)+
  #xlim(-2, 2.3)+
  xlab("PC1 (48.9%)")+
  ylab("PC2 (14.5%)")

# Water Canyon
trait_WC <- trait_master %>%
  filter(Population == "Water Canyon") %>%
  as.data.frame(decostand(., "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
trait_WC_matrix <- as.matrix(trait_WC[,7:ncol(trait_WC)])
pca_WC_trait = rda(trait_WC_matrix , scale = TRUE) #run PCA on all traits
biplot(pca_WC_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores_WC <- as.data.frame(scores(pca_WC_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab_WC = as.data.frame(cbind(trait_WC[,1:6],pca_trait_scores_WC)) #add plot info back
pca_trait_scores_lab_WC$Treatment <- apply(pca_trait_scores_lab_WC[ ,3:4 ] , 1 , paste , collapse = "_" )
envout<-as.data.frame(scores(pca_WC_trait , choices=c(1,2), display=c("species")))
summary(pca_WC_trait)
PCA_WC <- ggplot(pca_trait_scores_lab_WC, aes(x = PC1, y = PC2))+
  geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15))+
  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
            fontface = "bold", label = row.names(envout), size = 5)+
  #xlim(-2, 2.3)+
  xlab("PC1 (41.2%)")+
  ylab("PC2 (20.7%)")

# Reno
trait_Reno <- trait_master %>%
  filter(Population == "Reno") %>%
  as.data.frame(decostand(., "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
trait_Reno_matrix <- as.matrix(trait_Reno[,7:ncol(trait_Reno)])
pca_Reno_trait = rda(trait_Reno_matrix , scale = TRUE) #run PCA on all traits
biplot(pca_Reno_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores_Reno <- as.data.frame(scores(pca_Reno_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab_Reno = as.data.frame(cbind(trait_Reno[,1:6],pca_trait_scores_Reno)) #add plot info back
pca_trait_scores_lab_Reno$Treatment <- apply(pca_trait_scores_lab_Reno[ ,3:4 ] , 1 , paste , collapse = "_" )
envout<-as.data.frame(scores(pca_Reno_trait , choices=c(1,2), display=c("species")))
summary(pca_Reno_trait)
PCA_Reno <- ggplot(pca_trait_scores_lab_Reno, aes(x = PC1, y = PC2))+
  geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15))+
  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
            fontface = "bold", label = row.names(envout), size = 5)+
  #xlim(-2, 2.3)+
  xlab("PC1 (54.3%)")+
  ylab("PC2 (17.9%)")

ggarrange(PCA_Butte, PCA_Steens, PCA_EOARC, PCA_WC, PCA_Reno,
          common.legend = TRUE,
          labels = c("Butte Valley", "Steens", "EOARC", "Water Canyon", "Reno"), label.x = .1, label.y = .99, font.label = c(color = "blue"))

# PCA within each treatment:
# BRTE_DRY
trait_BRTE_dry <- trait_master %>%
  filter(Competition == "BRTE" & Water == "Dry") %>%
  as.data.frame(decostand(., "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
trait_BRTE_dry_matrix <- as.matrix(trait_BRTE_dry[,7:ncol(trait_BRTE_dry)])
pca_BRTE_dry_trait = rda(trait_BRTE_dry_matrix , scale = TRUE) #run PCA on all traits
biplot(pca_BRTE_dry_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores_BRTE_dry <- as.data.frame(scores(pca_BRTE_dry_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab_BRTE_dry = as.data.frame(cbind(trait_BRTE_dry[,1:6],pca_trait_scores_BRTE_dry)) #add plot info back
pca_trait_scores_lab_BRTE_dry$Population <- ordered(as.factor(pca_trait_scores_lab_BRTE_dry$Population), levels = c("Butte Valley","Steens","EOARC",
                                                                                                  "Water Canyon",  "Reno"))
envout<-as.data.frame(scores(pca_BRTE_dry_trait , choices=c(1,2), display=c("species")))
summary(pca_BRTE_dry_trait)
PCA_BRTE_dry <- ggplot(pca_trait_scores_lab_BRTE_dry, aes(x = PC1, y = PC2))+
                  geom_point(size = 4, aes(colour = Population), alpha = 0.5)+
                  theme(text = element_text(size=18),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                        axis.title = element_text(size = 15))+
                  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
                  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                               alpha = 0.5, size = 1, colour = "grey30") +
                  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
                            fontface = "bold", label = row.names(envout), size = 5)+
                  #xlim(-2, 2.3)+
                  xlab("PC1 (35.4%)")+
                  ylab("PC2 (17.9%)")

#BRTE_WET
trait_BRTE_wet <- trait_master %>%
  filter(Competition == "BRTE" & Water == "Wet") %>%
  as.data.frame(decostand(., "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
trait_BRTE_wet_matrix <- as.matrix(trait_BRTE_wet[,7:ncol(trait_BRTE_wet)])
pca_BRTE_wet_trait = rda(trait_BRTE_wet_matrix , scale = TRUE) #run PCA on all traits
biplot(pca_BRTE_wet_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores_BRTE_wet <- as.data.frame(scores(pca_BRTE_wet_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab_BRTE_wet = as.data.frame(cbind(trait_BRTE_wet[,1:6],pca_trait_scores_BRTE_wet)) #add plot info back
pca_trait_scores_lab_BRTE_wet$Population <- ordered(as.factor(pca_trait_scores_lab_BRTE_wet$Population), levels = c("Butte Valley","Steens","EOARC",
                                                                                                                    "Water Canyon",  "Reno"))
envout<-as.data.frame(scores(pca_BRTE_wet_trait , choices=c(1,2), display=c("species")))
summary(pca_BRTE_wet_trait)
PCA_BRTE_wet <- ggplot(pca_trait_scores_lab_BRTE_wet, aes(x = PC1, y = PC2))+
                  geom_point(size = 4, aes(colour = Population), alpha = 0.5)+
                  theme(text = element_text(size=18),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                        axis.title = element_text(size = 15))+
                  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
                  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                               alpha = 0.5, size = 1, colour = "grey30") +
                  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
                            fontface = "bold", label = row.names(envout), size = 5)+
                  #xlim(-2, 2.3)+
                  xlab("PC1 (31.1%)")+
                  ylab("PC2 (22.3%)")

# None_DRY
trait_None_dry <- trait_master %>%
  filter(Competition == "None" & Water == "Dry") %>%
  as.data.frame(decostand(., "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
trait_None_dry_matrix <- as.matrix(trait_None_dry[,7:ncol(trait_None_dry)])
pca_None_dry_trait = rda(trait_None_dry_matrix , scale = TRUE) #run PCA on all traits
biplot(pca_None_dry_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores_None_dry <- as.data.frame(scores(pca_None_dry_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab_None_dry = as.data.frame(cbind(trait_None_dry[,1:6],pca_trait_scores_None_dry)) #add plot info back
pca_trait_scores_lab_None_dry$Population <- ordered(as.factor(pca_trait_scores_lab_None_dry$Population), levels = c("Butte Valley","Steens","EOARC",
                                                                                                                    "Water Canyon",  "Reno"))
envout<-as.data.frame(scores(pca_None_dry_trait , choices=c(1,2), display=c("species")))
summary(pca_None_dry_trait)
PCA_None_dry <- ggplot(pca_trait_scores_lab_None_dry, aes(x = PC1, y = PC2))+
                  geom_point(size = 4, aes(colour = Population), alpha = 0.5)+
                  theme(text = element_text(size=18),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                        axis.title = element_text(size = 15))+
                  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
                  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                               alpha = 0.5, size = 1, colour = "grey30") +
                  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
                            fontface = "bold", label = row.names(envout), size = 5)+
                  #xlim(-2, 2.3)+
                  xlab("PC1 (29.3%)")+
                  ylab("PC2 (21.2%)")


#NONE_WET
trait_None_wet <- trait_master %>%
  filter(Competition == "None" & Water == "Wet") %>%
  as.data.frame(decostand(., "standardize")) %>%
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
trait_None_wet_matrix <- as.matrix(trait_None_wet[,7:ncol(trait_None_wet)])
pca_None_wet_trait = rda(trait_None_wet_matrix , scale = TRUE) #run PCA on all traits
biplot(pca_None_wet_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores_None_wet <- as.data.frame(scores(pca_None_wet_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab_None_wet = as.data.frame(cbind(trait_None_wet[,1:6],pca_trait_scores_None_wet)) #add plot info back
pca_trait_scores_lab_None_wet$Population <- ordered(as.factor(pca_trait_scores_lab_None_wet$Population), levels = c("Butte Valley","Steens","EOARC",
                                                                                                                    "Water Canyon",  "Reno"))
envout<-as.data.frame(scores(pca_None_wet_trait , choices=c(1,2), display=c("species")))
summary(pca_None_wet_trait)
PCA_None_wet <- ggplot(pca_trait_scores_lab_None_wet, aes(x = PC1, y = PC2))+
                  geom_point(size = 4, aes(colour = Population), alpha = 0.5)+
                  theme(text = element_text(size=18),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                        axis.title = element_text(size = 15))+
                  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
                  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                               alpha = 0.5, size = 1, colour = "grey30") +
                  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
                            fontface = "bold", label = row.names(envout), size = 5)+
                  #xlim(-2, 2.3)+
                  xlab("PC1 (36.9%)")+
                  ylab("PC2 (21.6%)")
ggarrange(PCA_None_wet, PCA_None_dry, PCA_BRTE_wet, PCA_BRTE_dry, 
          common.legend = TRUE,
          labels = c("None-Wet", "None-Dry", "BRTE-Wet", "BRTE-Dry"), label.x = .05, label.y = .99, font.label = c(color = "blue"))
#----------------------------------------------------------#
# How did each trait respond to each treatment?
# Calculate mean and standard error of each trait by population and treatment
trait_long <- cbind(trait_master[,1:6],decostand(trait_master[,7:ncol(trait_master)], "standardize")) %>%
  pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value")
trait_long$Population <- ordered(as.factor(trait_long$Population), levels = c("Butte Valley","Steens","EOARC",
                                                                              "Water Canyon",  "Reno"))
trait_long$Treatment <- apply(trait_long[ ,3:4 ] , 1 , paste , collapse = "-" )
trait_long$Treatment <- ordered(as.factor(trait_long$Treatment), levels = c("None-Wet", "None-Dry", "BRTE-Wet", "BRTE-Dry"))

mean_trait_long <- trait_long%>%
  group_by(Population, Treatment, trait) %>%
  summarise(mean = mean(value),
            se = se(value))
f_AG_trait <-ggplot(mean_trait_long%>%filter(trait%in%c("Emergence", "Height", "SLA", "LDMC")), aes(x = Population, y = mean, col = Treatment)) +
                geom_point(position = position_dodge(width = 0.5))+
                geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                facet_wrap(~trait, scales = "free", ncol = 2)+
                theme_bw()+
                scale_color_manual(values=c("#56B4E9","#E69F00", "#6A0DAD", "#999999" ))+
                ylab("")
f_BG_trait <-ggplot(mean_trait_long%>%filter(trait%in%c("RMR",  "Tips", "PropF", "SRL")), aes(x = Population, y = mean, col = Treatment)) +
                geom_point(position = position_dodge(width = 0.5))+
                geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                facet_wrap(~trait, scales = "free", ncol = 2)+
                theme_bw()+
                scale_color_manual(values=c("#56B4E9","#E69F00", "#6A0DAD", "#999999" ))+
                ylab("")

# Graph them together
ggarrange(f_AG_trait, f_BG_trait, ncol = 1, nrow = 2, labels = c("(a)", "(b)"),
          font.label = list(size = 15), common.legend = TRUE, legend = "right", heights = c(1, 1))
#----------------------------------------------------------#
# Which trait is correlated with higher POSE performance?
trait_standard <- cbind(trait_master[,1:6],decostand(trait_master[,7:ncol(trait_master)], "standardize"))

survival_trait <- inner_join(trait_standard, stemcount) %>%
  select(-POSE_emergence_stem_count) %>%
  select(-BRTE_stem_count) %>%
  drop_na()
survival_trait$Treatment <- apply(survival_trait[ ,3:4 ] , 1 , paste , collapse = "_" )

ggplot(survival_trait%>%select(-TotalBiomass)%>%
         pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value")%>%
         filter(trait%in%c("Emergence", "Height", "SLA", "LDMC", "RMR",  "Tips", "PropF", "SRL")), 
       aes(x = value, y = POSE_survival_stem_count/25, col = as.factor(Treatment)))+
  geom_jitter()+
  facet_wrap(~trait, scale = "free", ncol = 4)+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 12))+
  ylab(bquote(italic(P.~secunda)~Establishment~Rate)) +
  xlab("Trait z-scores")

totalbiomass_trait <- trait_standard %>%
  select(-TotalBiomass) %>%
  cbind(., trait_master[,16])
totalbiomass_trait$Treatment <- apply(totalbiomass_trait[ ,3:4 ] , 1 , paste , collapse = "_" )

ggplot(totalbiomass_trait %>%
  pivot_longer(cols = Length:LDMC, names_to = "trait", values_to = "value") %>%
  filter(trait%in%c("Emergence", "Height", "SLA", "LDMC", "RMR",  "Tips", "PropF", "SRL")),
  aes(x = value, y = TotalBiomass, col = Treatment))+
  geom_jitter()+
  facet_wrap(~trait, scale = "free", ncol = 4)+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 12))+
  ylab(bquote(italic(P.~secunda)~Total~Biomass~(g))) +
  xlab("Trait z-scores")

#----------------------------------------------------------#
# Dispersion and centroid of traits
# multivariate_difference (codyn) uses a Bray-Curtis dissimilarity matrix to calculate changes in composition (distance between centroids) and dispersion (distance around centroids)
trait_long_select <- trait_long %>%
  filter(trait%in%c("Emergence", "Height", "SLA", "LDMC", "RMR",  "Tips", "PropF", "SRL")) %>%
  mutate(value = value+3) #make all z-scores positive by adding 3

multivariate_difference(trait_long_select, 
                    time.var = NULL,
                    replicate.var = "PotID",
                    treatment.var = "Treatment",
                    species.var = "trait",
                    abundance.var = "value") #Overall, None-Wet and BRTE-Dry has the largest difference in trait composition and dispersion
multivariate_difference(trait_long_select, 
                        time.var = "Population",
                        replicate.var = "PotID",
                        treatment.var = "Treatment",
                        species.var = "trait",
                        abundance.var = "value",
                        reference.treatment = "None-Wet") #Reno has the largest difference in trait composition between None-Wet and BRTE-DRY.
adonis(value~Treatment, data = trait_long_select, perm = 999, method = "euclidean") #PERMANOVA results: Significant treatment effect p = 0.001

treatment <- apply(trait_master[ ,3:4 ] , 1 , paste , collapse = "_" )
dis <- vegdist(trait_matrix, method = "euclidean")
mod <- betadisper(dis, group = treatment)
anova(mod)
TukeyHSD(mod)
view(mod$vectors)

#----------------------------------------------------------#
# Is trait plasticity correlated with drought tolerance?

# Calculate PCA score distance as a metric of plasticity
plasticity <- pca_trait_scores_lab %>%
  dplyr::select(Population, Replicate, Competition, Water, PC1, PC2) %>%
  pivot_wider(names_from = Water, values_from = c(PC1, PC2)) %>%
  group_by(Population, Replicate, Competition) %>%
  na.omit()%>%
  summarise(d_trait = sqrt((PC1_Dry-PC1_Wet)^2 + (PC2_Dry- PC2_Wet)^2))

# Calculate drought tolerance: relative change in survival rate and biomass 
survival_delta <- growth%>%
  dplyr::select(Population, Replicate, Water, Competition, POSE_survival_stem_count) %>%
  mutate(survivalrate = POSE_survival_stem_count/25) %>%
  dplyr::select(Population, Replicate, Water, Competition, survivalrate) %>%
  pivot_wider(names_from = Water, values_from = survivalrate) %>%
  group_by(Population, Replicate, Competition) %>%
  summarise(survival_delta = (Dry-Wet)/Wet)
biomass_delta <- biomass_dat %>%
  filter(Life_stage == "seedling") %>%
  dplyr::select(Population, Replicate, Water, Competition, TotalBiomass) %>%
  pivot_wider(names_from = Water, values_from = TotalBiomass) %>%
  group_by(Population, Replicate, Competition) %>%
  na.omit()%>%
  summarise(biomass_delta = (Dry-Wet)/Wet)

# Combine dataframes to plot it
plastic_combined <- left_join(survival_delta, biomass_delta) %>%
  left_join(., plasticity) 
f1 <- ggplot(plastic_combined, aes(x = d_trait, y = survival_delta, col = Competition)) +
  geom_point(aes(col = Competition))+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 12))+
  ylab(bquote(Relative~Change~Establishment~Rate)) +
  xlab("Trait plasticity")+
  geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5)
f2 <- ggplot(plastic_combined, aes(x = d_trait, y = biomass_delta, col = Competition)) +
  geom_point(aes(col = Competition))+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 12))+
  ylab(bquote(Relative~Change~Total~Biomass)) +
  xlab("Trait plasticity")+
  geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5)


# Graph them together
ggarrange(f1, f2, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
          font.label = list(size = 15), common.legend = TRUE, legend = "right", label.x = -.02, label.y = .99)

# Alternative:
# Calculate trait plasticity or deviation from None_Wet as control 
plasticity_v2 <- pca_trait_scores_lab %>%
  dplyr::select(Population, Replicate, Treatment, PC1, PC2) %>%
  pivot_wider(names_from = Treatment, values_from = c(PC1, PC2)) %>%
  group_by(Population, Replicate) %>%
  na.omit()%>%
  summarise(None_Dry = sqrt((PC1_None_Dry-PC1_None_Wet)^2 + (PC2_None_Dry- PC2_None_Wet)^2),
            BRTE_Wet = sqrt((PC1_BRTE_Wet-PC1_None_Wet)^2 + (PC2_BRTE_Wet- PC2_None_Wet)^2),
            BRTE_Dry = sqrt((PC1_BRTE_Dry-PC1_None_Wet)^2 + (PC2_BRTE_Dry- PC2_None_Wet)^2)) %>%
  pivot_longer(cols = None_Dry:BRTE_Dry, names_to = "Treatment", values_to = "d_trait")

ggplot(plasticity_v2, aes(x = Population, y = d_trait, fill = Treatment)) +
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values=c("#999999", "#6A0DAD", "#E69F00"))

# Calculate stress tolerance: relative change in survival rate and biomass 
growth$Treatment <- apply(growth[ ,3:4 ] , 1 , paste , collapse = "_" )
biomass_dat$Treatment <- apply(biomass_dat[ ,3:4 ] , 1 , paste , collapse = "_" )
survival_delta_v2 <- growth%>%
  dplyr::select(Population, Replicate, Treatment, POSE_survival_stem_count) %>%
  mutate(survivalrate = POSE_survival_stem_count/25) %>%
  dplyr::select(Population, Replicate, Treatment, survivalrate) %>%
  pivot_wider(names_from = Treatment, values_from = survivalrate) %>%
  group_by(Population, Replicate) %>%
  summarise(None_Dry = (None_Dry-None_Wet)/None_Wet,
            BRTE_Wet = (BRTE_Wet-None_Wet)/None_Wet,
            BRTE_Dry = (BRTE_Dry-None_Wet)/None_Wet) %>%
  pivot_longer(cols = None_Dry:BRTE_Dry, names_to = "Treatment", values_to = "survival_delta")
biomass_delta_v2 <- biomass_dat %>%
  filter(Life_stage == "seedling") %>%
  dplyr::select(Population, Replicate, Treatment, TotalBiomass) %>%
  pivot_wider(names_from = Treatment, values_from = TotalBiomass) %>%
  group_by(Population, Replicate) %>%
  na.omit()%>%
  summarise(None_Dry = (None_Dry-None_Wet)/None_Wet,
            BRTE_Wet = (BRTE_Wet-None_Wet)/None_Wet,
            BRTE_Dry = (BRTE_Dry-None_Wet)/None_Wet) %>%
  pivot_longer(cols = None_Dry:BRTE_Dry, names_to = "Treatment", values_to = "biomass_delta")

# Combine dataframes to plot it
plastic_combined_v2 <- left_join(survival_delta_v2, biomass_delta_v2) %>%
  left_join(., plasticity_v2) 
f1_V2 <- ggplot(plastic_combined_v2, aes(x = d_trait, y = survival_delta)) +
  geom_point(aes(col = Treatment))+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 12))+
  ylab(bquote(Relative~Change~Establishment~Rate)) +
  xlab("Trait plasticity")+
  geom_smooth(method = "lm", formula = y~poly(x,3), size=0.5, color = "black")+
  scale_color_manual(values=c("#999999", "#6A0DAD", "#E69F00"))
f2_v2 <- ggplot(plastic_combined_v2, aes(x = d_trait, y = biomass_delta)) +
  geom_point(aes(col = Treatment))+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 12))+
  ylab(bquote(Relative~Change~Total~Biomass)) +
  xlab("Trait plasticity")+
  geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5, color = "black")+
  scale_color_manual(values=c("#999999", "#6A0DAD", "#E69F00"))

# Graph them together
ggarrange(f1_V2, f2_v2, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
          font.label = list(size = 15), common.legend = TRUE, legend = "right", label.x = -.02, label.y = .99)

