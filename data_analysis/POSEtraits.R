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
library(ggrepel)

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
  dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -AvgDiam, -TotalBiomass)

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
  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1*0.9, yend = PC2*0.9),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text_repel(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
            fontface = "bold", label = row.names(envout), size = 5, force = 0.02)+
  xlim(-2, 2)+
  ylim(-1.6, 1.6)+
  xlab("PC1 (43.0%)")+
  ylab("PC2 (14.4%)")+
  scale_color_manual(values=c("#FBD947", "#FF8C07",  "#D6AACE","#CE026E"))

#calculate centroids within each treatment
pca.centroids <- pca_trait_scores_lab %>%
  group_by(Treatment) %>%
  summarise(mean_PC1 = median(PC1),
            mean_PC2 = median(PC2))

ggplot(pca.centroids, aes(x = mean_PC1, y = mean_PC2))+
  geom_point(size = 4, aes(col = Treatment), alpha = 0.5)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15))+
  xlim(-0.7, 0.7)+
  ylim(-0.7, 0.7)+
  xlab("PC1 (43.0%)")+
  ylab("PC2 (14.4%)")+
  scale_color_manual(values=c("#FBD947", "#FF9507", "#D6AACE","#CE026E"))

# calculate dispersion within each treatment
pca.dispersion <- pca_trait_scores_lab %>%
  group_by(Treatment) %>%
  summarise(dispersion = sqrt((PC1-median(PC1))^2 + (PC2-median(PC2))^2))
pca.dispersion.summary <- pca.dispersion %>%
  group_by(Treatment) %>%
  summarise(mean = mean(dispersion),
            se = se(dispersion))
ggplot(pca.dispersion.summary, aes(x = Treatment, y = mean, fill = Treatment))+
  geom_bar(stat = "identity" , show.legend = FALSE)+
  #facet_wrap(vars(Population), ncol = 5)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
  xlab("Treatment")+
  ylab("Trait Dispersion")+
  scale_fill_manual(values=c("#FBD947", "#FF9507",  "#D6AACE","#CE026E"))+
  scale_x_discrete(labels=c("BRTE_Dry" = "BD", "BRTE_Wet" = "BW",
                            "None_Dry" = "ND", "None_Wet" = "NW"))

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

# ---------------------------------------------------------------#  
#  PCA within each population:
# # Butte Valley
# trait_Butte <- trait_master %>%
#   filter(Population == "Butte Valley") %>%
#   as.data.frame(decostand(., "standardize")) %>%
#   dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
# trait_Butte_matrix <- as.matrix(trait_Butte[,7:ncol(trait_Butte)])
# pca_Butte_trait = rda(trait_Butte_matrix , scale = TRUE) #run PCA on all traits
# biplot(pca_Butte_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
# pca_trait_scores_Butte <- as.data.frame(scores(pca_Butte_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
# pca_trait_scores_lab_Butte = as.data.frame(cbind(trait_Butte[,1:6],pca_trait_scores_Butte)) #add plot info back
# pca_trait_scores_lab_Butte$Treatment <- apply(pca_trait_scores_lab_Butte[ ,3:4 ] , 1 , paste , collapse = "_" )
# envout<-as.data.frame(scores(pca_Butte_trait , choices=c(1,2), display=c("species")))
# summary(pca_Butte_trait)
# PCA_Butte <- ggplot(pca_trait_scores_lab_Butte, aes(x = PC1, y = PC2))+
#   geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
#   theme(text = element_text(size=18),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         axis.title = element_text(size = 15))+
#   scale_shape_manual(values=c(8, 5, 15, 17, 19))+
#   geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                alpha = 0.5, size = 1, colour = "grey30") +
#   geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
#             fontface = "bold", label = row.names(envout), size = 5)+
#   #xlim(-2, 2.3)+
#   xlab("PC1 (32.7%)")+
#   ylab("PC2 (21.9%)")
# 
# # Steens
# trait_Steens <- trait_master %>%
#   filter(Population == "Steens") %>%
#   as.data.frame(decostand(., "standardize")) %>%
#   dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
# trait_Steens_matrix <- as.matrix(trait_Steens[,7:ncol(trait_Steens)])
# pca_Steens_trait = rda(trait_Steens_matrix , scale = TRUE) #run PCA on all traits
# biplot(pca_Steens_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
# pca_trait_scores_Steens <- as.data.frame(scores(pca_Steens_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
# pca_trait_scores_lab_Steens = as.data.frame(cbind(trait_Steens[,1:6],pca_trait_scores_Steens)) #add plot info back
# pca_trait_scores_lab_Steens$Treatment <- apply(pca_trait_scores_lab_Steens[ ,3:4 ] , 1 , paste , collapse = "_" )
# envout<-as.data.frame(scores(pca_Steens_trait , choices=c(1,2), display=c("species")))
# summary(pca_Steens_trait)
# PCA_Steens <- ggplot(pca_trait_scores_lab_Steens, aes(x = PC1, y = PC2))+
#   geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
#   theme(text = element_text(size=18),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         axis.title = element_text(size = 15))+
#   scale_shape_manual(values=c(8, 5, 15, 17, 19))+
#   geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                alpha = 0.5, size = 1, colour = "grey30") +
#   geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
#             fontface = "bold", label = row.names(envout), size = 5)+
#   #xlim(-2, 2.3)+
#   xlab("PC1 (42.1%)")+
#   ylab("PC2 (18.4%)")
# 
# # EOARC
# trait_EOARC <- trait_master %>%
#   filter(Population == "EOARC") %>%
#   as.data.frame(decostand(., "standardize")) %>%
#   dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
# trait_EOARC_matrix <- as.matrix(trait_EOARC[,7:ncol(trait_EOARC)])
# pca_EOARC_trait = rda(trait_EOARC_matrix , scale = TRUE) #run PCA on all traits
# biplot(pca_EOARC_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
# pca_trait_scores_EOARC <- as.data.frame(scores(pca_EOARC_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
# pca_trait_scores_lab_EOARC = as.data.frame(cbind(trait_EOARC[,1:6],pca_trait_scores_EOARC)) #add plot info back
# pca_trait_scores_lab_EOARC$Treatment <- apply(pca_trait_scores_lab_EOARC[ ,3:4 ] , 1 , paste , collapse = "_" )
# envout<-as.data.frame(scores(pca_EOARC_trait , choices=c(1,2), display=c("species")))
# summary(pca_EOARC_trait)
# PCA_EOARC <- ggplot(pca_trait_scores_lab_EOARC, aes(x = PC1, y = PC2))+
#   geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
#   theme(text = element_text(size=18),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         axis.title = element_text(size = 15))+
#   scale_shape_manual(values=c(8, 5, 15, 17, 19))+
#   geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                alpha = 0.5, size = 1, colour = "grey30") +
#   geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
#             fontface = "bold", label = row.names(envout), size = 5)+
#   #xlim(-2, 2.3)+
#   xlab("PC1 (48.9%)")+
#   ylab("PC2 (14.5%)")
# 
# # Water Canyon
# trait_WC <- trait_master %>%
#   filter(Population == "Water Canyon") %>%
#   as.data.frame(decostand(., "standardize")) %>%
#   dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
# trait_WC_matrix <- as.matrix(trait_WC[,7:ncol(trait_WC)])
# pca_WC_trait = rda(trait_WC_matrix , scale = TRUE) #run PCA on all traits
# biplot(pca_WC_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
# pca_trait_scores_WC <- as.data.frame(scores(pca_WC_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
# pca_trait_scores_lab_WC = as.data.frame(cbind(trait_WC[,1:6],pca_trait_scores_WC)) #add plot info back
# pca_trait_scores_lab_WC$Treatment <- apply(pca_trait_scores_lab_WC[ ,3:4 ] , 1 , paste , collapse = "_" )
# envout<-as.data.frame(scores(pca_WC_trait , choices=c(1,2), display=c("species")))
# summary(pca_WC_trait)
# PCA_WC <- ggplot(pca_trait_scores_lab_WC, aes(x = PC1, y = PC2))+
#   geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
#   theme(text = element_text(size=18),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         axis.title = element_text(size = 15))+
#   scale_shape_manual(values=c(8, 5, 15, 17, 19))+
#   geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                alpha = 0.5, size = 1, colour = "grey30") +
#   geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
#             fontface = "bold", label = row.names(envout), size = 5)+
#   #xlim(-2, 2.3)+
#   xlab("PC1 (41.2%)")+
#   ylab("PC2 (20.7%)")
# 
# # Reno
# trait_Reno <- trait_master %>%
#   filter(Population == "Reno") %>%
#   as.data.frame(decostand(., "standardize")) %>%
#   dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
# trait_Reno_matrix <- as.matrix(trait_Reno[,7:ncol(trait_Reno)])
# pca_Reno_trait = rda(trait_Reno_matrix , scale = TRUE) #run PCA on all traits
# biplot(pca_Reno_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
# pca_trait_scores_Reno <- as.data.frame(scores(pca_Reno_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
# pca_trait_scores_lab_Reno = as.data.frame(cbind(trait_Reno[,1:6],pca_trait_scores_Reno)) #add plot info back
# pca_trait_scores_lab_Reno$Treatment <- apply(pca_trait_scores_lab_Reno[ ,3:4 ] , 1 , paste , collapse = "_" )
# envout<-as.data.frame(scores(pca_Reno_trait , choices=c(1,2), display=c("species")))
# summary(pca_Reno_trait)
# PCA_Reno <- ggplot(pca_trait_scores_lab_Reno, aes(x = PC1, y = PC2))+
#   geom_point(size = 4, aes(colour = Treatment), alpha = 0.5)+
#   theme(text = element_text(size=18),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         axis.title = element_text(size = 15))+
#   scale_shape_manual(values=c(8, 5, 15, 17, 19))+
#   geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                alpha = 0.5, size = 1, colour = "grey30") +
#   geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
#             fontface = "bold", label = row.names(envout), size = 5)+
#   #xlim(-2, 2.3)+
#   xlab("PC1 (54.3%)")+
#   ylab("PC2 (17.9%)")
# 
# ggarrange(PCA_Butte, PCA_Steens, PCA_EOARC, PCA_WC, PCA_Reno,
#           common.legend = TRUE,
#           labels = c("Butte Valley", "Steens", "EOARC", "Water Canyon", "Reno"), label.x = .1, label.y = .99, font.label = c(color = "blue"))
# 
# # PCA within each treatment:
# # BRTE_DRY
# trait_BRTE_dry <- trait_master %>%
#   filter(Competition == "BRTE" & Water == "Dry") %>%
#   as.data.frame(decostand(., "standardize")) %>%
#   dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
# trait_BRTE_dry_matrix <- as.matrix(trait_BRTE_dry[,7:ncol(trait_BRTE_dry)])
# pca_BRTE_dry_trait = rda(trait_BRTE_dry_matrix , scale = TRUE) #run PCA on all traits
# biplot(pca_BRTE_dry_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
# pca_trait_scores_BRTE_dry <- as.data.frame(scores(pca_BRTE_dry_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
# pca_trait_scores_lab_BRTE_dry = as.data.frame(cbind(trait_BRTE_dry[,1:6],pca_trait_scores_BRTE_dry)) #add plot info back
# pca_trait_scores_lab_BRTE_dry$Population <- ordered(as.factor(pca_trait_scores_lab_BRTE_dry$Population), levels = c("Butte Valley","Steens","EOARC",
#                                                                                                   "Water Canyon",  "Reno"))
# envout<-as.data.frame(scores(pca_BRTE_dry_trait , choices=c(1,2), display=c("species")))
# summary(pca_BRTE_dry_trait)
# PCA_BRTE_dry <- ggplot(pca_trait_scores_lab_BRTE_dry, aes(x = PC1, y = PC2))+
#                   geom_point(size = 4, aes(colour = Population), alpha = 0.5)+
#                   theme(text = element_text(size=18),
#                         panel.grid.major = element_blank(),
#                         panel.grid.minor = element_blank(),
#                         panel.background = element_blank(),
#                         axis.line = element_line(colour = "black"),
#                         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#                         axis.title = element_text(size = 15))+
#                   scale_shape_manual(values=c(8, 5, 15, 17, 19))+
#                   geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                                alpha = 0.5, size = 1, colour = "grey30") +
#                   geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
#                             fontface = "bold", label = row.names(envout), size = 5)+
#                   #xlim(-2, 2.3)+
#                   xlab("PC1 (35.4%)")+
#                   ylab("PC2 (17.9%)")
# 
# #BRTE_WET
# trait_BRTE_wet <- trait_master %>%
#   filter(Competition == "BRTE" & Water == "Wet") %>%
#   as.data.frame(decostand(., "standardize")) %>%
#   dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
# trait_BRTE_wet_matrix <- as.matrix(trait_BRTE_wet[,7:ncol(trait_BRTE_wet)])
# pca_BRTE_wet_trait = rda(trait_BRTE_wet_matrix , scale = TRUE) #run PCA on all traits
# biplot(pca_BRTE_wet_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
# pca_trait_scores_BRTE_wet <- as.data.frame(scores(pca_BRTE_wet_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
# pca_trait_scores_lab_BRTE_wet = as.data.frame(cbind(trait_BRTE_wet[,1:6],pca_trait_scores_BRTE_wet)) #add plot info back
# pca_trait_scores_lab_BRTE_wet$Population <- ordered(as.factor(pca_trait_scores_lab_BRTE_wet$Population), levels = c("Butte Valley","Steens","EOARC",
#                                                                                                                     "Water Canyon",  "Reno"))
# envout<-as.data.frame(scores(pca_BRTE_wet_trait , choices=c(1,2), display=c("species")))
# summary(pca_BRTE_wet_trait)
# PCA_BRTE_wet <- ggplot(pca_trait_scores_lab_BRTE_wet, aes(x = PC1, y = PC2))+
#                   geom_point(size = 4, aes(colour = Population), alpha = 0.5)+
#                   theme(text = element_text(size=18),
#                         panel.grid.major = element_blank(),
#                         panel.grid.minor = element_blank(),
#                         panel.background = element_blank(),
#                         axis.line = element_line(colour = "black"),
#                         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#                         axis.title = element_text(size = 15))+
#                   scale_shape_manual(values=c(8, 5, 15, 17, 19))+
#                   geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                                alpha = 0.5, size = 1, colour = "grey30") +
#                   geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
#                             fontface = "bold", label = row.names(envout), size = 5)+
#                   #xlim(-2, 2.3)+
#                   xlab("PC1 (31.1%)")+
#                   ylab("PC2 (22.3%)")
# 
# # None_DRY
# trait_None_dry <- trait_master %>%
#   filter(Competition == "None" & Water == "Dry") %>%
#   as.data.frame(decostand(., "standardize")) %>%
#   dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
# trait_None_dry_matrix <- as.matrix(trait_None_dry[,7:ncol(trait_None_dry)])
# pca_None_dry_trait = rda(trait_None_dry_matrix , scale = TRUE) #run PCA on all traits
# biplot(pca_None_dry_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
# pca_trait_scores_None_dry <- as.data.frame(scores(pca_None_dry_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
# pca_trait_scores_lab_None_dry = as.data.frame(cbind(trait_None_dry[,1:6],pca_trait_scores_None_dry)) #add plot info back
# pca_trait_scores_lab_None_dry$Population <- ordered(as.factor(pca_trait_scores_lab_None_dry$Population), levels = c("Butte Valley","Steens","EOARC",
#                                                                                                                     "Water Canyon",  "Reno"))
# envout<-as.data.frame(scores(pca_None_dry_trait , choices=c(1,2), display=c("species")))
# summary(pca_None_dry_trait)
# PCA_None_dry <- ggplot(pca_trait_scores_lab_None_dry, aes(x = PC1, y = PC2))+
#                   geom_point(size = 4, aes(colour = Population), alpha = 0.5)+
#                   theme(text = element_text(size=18),
#                         panel.grid.major = element_blank(),
#                         panel.grid.minor = element_blank(),
#                         panel.background = element_blank(),
#                         axis.line = element_line(colour = "black"),
#                         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#                         axis.title = element_text(size = 15))+
#                   scale_shape_manual(values=c(8, 5, 15, 17, 19))+
#                   geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                                alpha = 0.5, size = 1, colour = "grey30") +
#                   geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
#                             fontface = "bold", label = row.names(envout), size = 5)+
#                   #xlim(-2, 2.3)+
#                   xlab("PC1 (29.3%)")+
#                   ylab("PC2 (21.2%)")
# 
# 
# #NONE_WET
# trait_None_wet <- trait_master %>%
#   filter(Competition == "None" & Water == "Wet") %>%
#   as.data.frame(decostand(., "standardize")) %>%
#   dplyr::select( -Forks, -Coarse, -Fine, -SurfArea, -Length, -AvgDiam, -TotalBiomass)
# trait_None_wet_matrix <- as.matrix(trait_None_wet[,7:ncol(trait_None_wet)])
# pca_None_wet_trait = rda(trait_None_wet_matrix , scale = TRUE) #run PCA on all traits
# biplot(pca_None_wet_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
# pca_trait_scores_None_wet <- as.data.frame(scores(pca_None_wet_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
# pca_trait_scores_lab_None_wet = as.data.frame(cbind(trait_None_wet[,1:6],pca_trait_scores_None_wet)) #add plot info back
# pca_trait_scores_lab_None_wet$Population <- ordered(as.factor(pca_trait_scores_lab_None_wet$Population), levels = c("Butte Valley","Steens","EOARC",
#                                                                                                                     "Water Canyon",  "Reno"))
# envout<-as.data.frame(scores(pca_None_wet_trait , choices=c(1,2), display=c("species")))
# summary(pca_None_wet_trait)
# PCA_None_wet <- ggplot(pca_trait_scores_lab_None_wet, aes(x = PC1, y = PC2))+
#                   geom_point(size = 4, aes(colour = Population), alpha = 0.5)+
#                   theme(text = element_text(size=18),
#                         panel.grid.major = element_blank(),
#                         panel.grid.minor = element_blank(),
#                         panel.background = element_blank(),
#                         axis.line = element_line(colour = "black"),
#                         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#                         axis.title = element_text(size = 15))+
#                   scale_shape_manual(values=c(8, 5, 15, 17, 19))+
#                   geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                                alpha = 0.5, size = 1, colour = "grey30") +
#                   geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30",
#                             fontface = "bold", label = row.names(envout), size = 5)+
#                   #xlim(-2, 2.3)+
#                   xlab("PC1 (36.9%)")+
#                   ylab("PC2 (21.6%)")
# ggarrange(PCA_None_wet, PCA_None_dry, PCA_BRTE_wet, PCA_BRTE_dry, 
#           common.legend = TRUE,
#           labels = c("None-Wet", "None-Dry", "BRTE-Wet", "BRTE-Dry"), label.x = .05, label.y = .99, font.label = c(color = "blue"))
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
                scale_color_manual(values=c("#CE026E","#D6AACE","#FF9507","#FBD947" ))+
                ylab("")
f_BG_trait <-ggplot(mean_trait_long%>%filter(trait%in%c("RMR",  "Tips", "PropF", "SRL")), aes(x = Population, y = mean, col = Treatment)) +
                geom_point(position = position_dodge(width = 0.5))+
                geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                facet_wrap(~trait, scales = "free", ncol = 2)+
                theme_bw()+
                scale_color_manual(values=c("#CE026E","#D6AACE","#FF9507","#FBD947"))+
                ylab("")

# Graph them together
ggarrange(f_AG_trait, f_BG_trait, ncol = 1, nrow = 2, labels = c("(a)", "(b)"),
          font.label = list(size = 15), common.legend = TRUE, legend = "right", heights = c(1, 1))
#----------------------------------------------------------#
# Which trait is correlated with higher POSE performance?
trait_standard <- cbind(trait_master[,1:6],decostand(trait_master[,7:ncol(trait_master)], "standardize"))

survival_trait <- inner_join(trait_standard, stemcount) %>%
  dplyr::select(-POSE_emergence_stem_count) %>%
  #dplyr::select(-BRTE_stem_count) %>%
  drop_na()
survival_trait$Treatment <- apply(survival_trait[ ,3:4 ] , 1 , paste , collapse = "_" )
survival_trait$Population <- ordered(as.factor(survival_trait$Population), levels = c("Butte Valley","Steens","EOARC",
                                                                              "Water Canyon",  "Reno"))

summary(lm(POSE_survival_stem_count ~ Emergence, data = survival_trait%>%filter(Competition == "BRTE"))) #0.2052, 3.393e-05
summary(lm(POSE_survival_stem_count ~ Emergence, data = survival_trait%>%filter(Competition == "None"))) #0.2052, 3.393e-05
summary(lm(POSE_survival_stem_count ~ Height, data = survival_trait%>%filter(Competition == "BRTE"))) #NA
summary(lm(POSE_survival_stem_count ~ Height, data = survival_trait%>%filter(Competition == "None"))) #NA
summary(lm(POSE_survival_stem_count ~ LDMC, data = survival_trait%>%filter(Competition == "BRTE"))) #0.05798 , 0.02262
summary(lm(POSE_survival_stem_count ~ LDMC, data = survival_trait%>%filter(Competition == "None"))) #NA
summary(lm(POSE_survival_stem_count ~ Length, data = survival_trait%>%filter(Competition == "BRTE"))) #0.2079, 2.999e-05
summary(lm(POSE_survival_stem_count ~ Length, data = survival_trait%>%filter(Competition == "None"))) #0.2079, 1.437e-05
summary(lm(POSE_survival_stem_count ~ PropF, data = survival_trait%>%filter(Competition == "BRTE"))) #0.4314, 1.678e-10
summary(lm(POSE_survival_stem_count ~ PropF, data = survival_trait%>%filter(Competition == "None"))) #0.2114, 1.201e-05
summary(lm(POSE_survival_stem_count ~ RMR, data = survival_trait%>%filter(Competition == "BRTE"))) #0.0881, 0.006208
summary(lm(POSE_survival_stem_count ~ RMR, data = survival_trait%>%filter(Competition == "None"))) #0.05738, 0.01893
summary(lm(POSE_survival_stem_count ~ SLA, data = survival_trait%>%filter(Competition == "BRTE"))) #0.05724, 0.02335
summary(lm(POSE_survival_stem_count ~ SLA, data = survival_trait%>%filter(Competition == "None"))) #NA
summary(lm(POSE_survival_stem_count ~ SRL, data = survival_trait%>%filter(Competition == "BRTE"))) #0.2617, 2.236e-06
summary(lm(POSE_survival_stem_count ~ SRL, data = survival_trait%>%filter(Competition == "None"))) #0.2718, 5.031e-07
summary(lm(POSE_survival_stem_count ~ Tips, data = survival_trait%>%filter(Competition == "BRTE"))) #0.04892, 0.03345
summary(lm(POSE_survival_stem_count ~ Tips, data = survival_trait%>%filter(Competition == "None"))) #0.06483, 0.0134

ann_text_traits_1 <- data.frame(x = c(2.5, 1.5, 0.9, 0, 0, 3.2, 2, 0.5),
                               y = 0.8,
                               label = c("BRTE R2 = 0.20, p < 0.001", "BRTE R2 = 0.05, p = 0.02", "BRTE R2 = 0.21, p < 0.001",
                                         "BRTE R2 = 0.43, p < 0.001", "BRTE R2 = 0.08, p = 0.006", "BRTE R2 = 0.05, p = 0.02",
                                         "BRTE R2 = 0.26, p < 0.001", "BRTE R2 = 0.04, p = 0.03"),
                               trait =  c("Emergence", "LDMC", "Length", "PropF", "RMR", "SLA", "SRL", "Tips"))
ann_text_traits_2 <- data.frame(x = c(2.5, 0.9, 0, 0, 2, 0.5),
                                y = 0.7,
                                label = c("None R2 = 0.20, p < 0.001", "None R2 = 0.05, p < 0.001", "None R2 = 0.21, p < 0.001",
                                          "None R2 = 0.05, p = 0.01", "None R2 = 0.27, p < 0.001", "None R2 = 0.06, p = 0.01"),
                                trait =  c("Emergence", "Length", "PropF", "RMR", "SRL", "Tips"))

trait_names <- list(
  'Emergence' = "Emergence Date",
  'Height' = "Height",
  'LDMC' = "Leaf Dry Matter Content",
  'Length' = "Root Length",
  'PropF' = "Proportion Fine Roots",
  'RMR' = "Root Mass Ratio",
  'SLA' = "Specific Leaf Area",
  'SRL' = "Specific Root Length",
  'Tips' = "Root Tips"
)
trait_labeller <- function(variable,value){
  return(trait_names[value])
}
ggplot(survival_trait%>%dplyr::select(-TotalBiomass)%>%
         pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value")%>%
         filter(trait%in%c("Emergence", "Height", "SLA", "LDMC", "RMR",  "Tips", "PropF", "SRL", "Length")), 
       aes(x = value, y = POSE_survival_stem_count/25))+
  geom_jitter(aes(col = Treatment))+
  facet_wrap(~trait, scale = "free", ncol = 3,  labeller=trait_labeller)+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 12))+
  ylab(bquote(Establishment~Rate)) +
  xlab("Trait z-scores")+
  scale_color_manual(values=c("#FBD947", "#FF9507",  "#D6AACE","#CE026E"))+
  geom_smooth(data = survival_trait%>%dplyr::select(-TotalBiomass)%>%
                pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value")
              %>%filter(trait %in% c("Emergence","LDMC", "Length", "PropF", "RMR","SLA", "SRL", "Tips")), 
              aes(linetype = Competition), col = "black",method = "lm",  size=0.9, se = FALSE)+
  geom_text(data = ann_text_traits_1, mapping = aes(x = x, y = y, label = label))+
  geom_text(data = ann_text_traits_2, mapping = aes(x = x, y = y, label = label))#+
  # geom_smooth(data = survival_trait%>%dplyr::select(-TotalBiomass)%>%
  #               pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value")
  #             %>%filter(trait %in% c("LDMC", "SLA"), 
  #              col = "black",method = "lm",  size=0.9, se = FALSE))

# ggplot(survival_trait%>%dplyr::select(-TotalBiomass)%>%filter(Competition == "BRTE") %>%
#          pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value")%>%
#          filter(trait%in%c("Emergence", "Height", "SLA", "LDMC", "RMR",  "Tips", "PropF", "SRL", "Length")), 
#        aes(y = value, x = BRTE_stem_count))+
#   geom_jitter(aes(col = Population))+
#   facet_wrap(~trait, scale = "free", ncol = 3)+
#   theme(text = element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#         axis.title = element_text(size = 12))+
#   ylab(bquote(trait)) +
#   xlab("BRTE stem") #+
#   #geom_smooth(aes(linetype = Competition), col = "black",method = "lm",  size=0.5)+
#   #scale_color_manual(values=c("#FBD947", "#FF9507",  "#86BBE8","#0240FF"))
# 
ggplot(survival_trait%>%dplyr::select(-TotalBiomass)%>%
         pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value")%>%
         filter(trait%in%c("Tips")),
       aes(x = value, y = POSE_survival_stem_count/25))+
  geom_jitter(aes(col = Treatment))+
  facet_wrap(~trait+Population, scale = "free", ncol = 5)+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 12))+
  ylab(bquote(italic(P.~secunda)~Establishment~Rate)) +
  xlab("Trait z-scores")+
  geom_smooth(aes(linetype = Competition), col = "black",method = "lm",  size=0.5)+
  scale_color_manual(values=c("#FBD947", "#FF9507",  "#86BBE8","#0240FF"))

totalbiomass_trait <- trait_standard %>%
  select(-TotalBiomass) %>%
  cbind(., trait_master[,16])
totalbiomass_trait$Treatment <- apply(totalbiomass_trait[ ,3:4 ] , 1 , paste , collapse = "_" )

ggplot(totalbiomass_trait %>%
  pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value") %>%
  filter(trait%in%c("Emergence", "Height", "SLA", "LDMC", "RMR",  "Tips", "PropF", "SRL",  "Length"))%>%
    filter(Treatment%in%c("None_Dry", "None_Wet")),
  aes(x = value, y = TotalBiomass))+
  geom_jitter(aes(col = Treatment))+
  facet_wrap(~trait, scale = "free", ncol = 3)+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 12))+
  ylab(bquote(italic(P.~secunda)~Total~Biomass~(g))) +
  xlab("Trait z-scores")+
  geom_smooth(aes(linetype = Competition), col = "black",method = "lm",  size=0.5)+
  scale_color_manual(values=c("#FBD947", "#FDAB02",   "#D6AACE","#CE026E"))

#----------------------------------------------------------#
# Does trait shift hurt or help POSE survival?
# plot PC axes and survival rate
PCA.survival <- inner_join(pca_trait_scores_lab, stemcount)
PC1_survival <- ggplot(PCA.survival, aes(x = PC1, y = POSE_survival_stem_count,  col = Treatment))+
  geom_point()+
  facet_wrap(vars(Population), ncol = 1)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.title = element_blank())+
  scale_color_manual(values=c("#FBD947", "#FF9507",  "#86BBE8","#0240FF"))
PC2_survival <- ggplot(PCA.survival, aes(x = PC2, y = POSE_survival_stem_count,  col = Treatment))+
  geom_point()+
  facet_wrap(vars(Population), ncol = 1)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.title = element_blank())+
  scale_color_manual(values=c("#FBD947", "#FF9507",  "#86BBE8","#0240FF"))
ggarrange(PC1_survival, PC2_survival, common.legend = TRUE)

summary.PCA.survival <- PCA.survival %>%
  group_by(Population,  Water, Competition) %>%
  summarise(mean_PC1 = mean(PC1),
            mean_PC2 = mean(PC2),
            mean_POSE = mean(POSE_survival_stem_count))
delta.PCA.survival <- summary.PCA.survival %>%
  group_by(Population, Competition) %>%
  summarise(delta_PC1 = mean_PC1[Water =="Dry"]-mean_PC1[Water == "Wet"],
            delta_PC2 = mean_PC2[Water == "Dry"]-mean_PC2[Water == "Wet"],
            delta_POSE = (mean_POSE[Water == "Dry"]-mean_POSE[Water == "Wet"])/mean_POSE[Water == "Wet"],
            distance = sqrt((mean_PC1[Water == "Dry"]-mean_PC1[Water == "Wet"])^2 + (mean_PC2[Water == "Dry"]-mean_PC2[Water == "Wet"])^2))

ggplot(delta.PCA.survival, aes(x = distance, y = delta_POSE, col = Competition))+
                            geom_point( size = 3)+
                            theme(text = element_text(size=18),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.background = element_blank(),
                                  axis.line = element_line(colour = "black"),
                                  panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                                  axis.title = element_text(size = 15))+
                            geom_vline(xintercept = 0, linetype = "dashed")+
                            geom_hline(yintercept = 0, linetype = "dashed")+
                            geom_text(aes(label=Population),hjust=-0.1, vjust=1)+
                            ylim( -1, 0.1)+
                            xlim(0, 1.1)+
                            ylab(bquote(Delta~Establishment~Rate))+
                            xlab(bquote(Euclidean~Distance))+
                            scale_color_manual(values=c("#FF9507","#0240FF"))
distance_None <- ggplot(delta.PCA.survival %>% filter(Competition == "None"), aes(x = distance, y = delta_POSE))+
                      geom_point(colour = "#111111", size = 3)+
                      theme(text = element_text(size=18),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                            axis.title = element_text(size = 15),
                            legend.title = element_blank())+
                      geom_vline(xintercept = 0, linetype = "dashed")+
                      geom_hline(yintercept = 0, linetype = "dashed")+
                      geom_text(aes(label=Population),hjust=-0.1, vjust=1)+
                      ylim( -6.2, 1)+
                      xlim(0, 1.1)+
                      ylab(bquote(Delta~Establishment~Rate))+
                      xlab(bquote(Euclidean~Distance))

distance_BRTE <- ggplot(delta.PCA.survival %>% filter(Competition == "BRTE"), aes(x = distance, y = delta_POSE))+
                      geom_point(colour = "#777777", size = 3)+
                      theme(text = element_text(size=18),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                            axis.title = element_text(size = 15),
                            legend.title = element_blank())+
                      geom_vline(xintercept = 0, linetype = "dashed")+
                      geom_hline(yintercept = 0, linetype = "dashed")+
                      geom_text(aes(label=Population),hjust=-0.1, vjust=1)+
                      ylim( -6.2, 1)+
                      xlim(-0, 1.1)+
                      ylab(bquote(Delta~Establishment~Rate))+
                      xlab(bquote(Euclidean~Distance))
delta_PC1_None <- ggplot(delta.PCA.survival%>% filter(Competition == "None"), aes(x = delta_PC1, y = delta_POSE))+
                      geom_point(colour = "#111111", size = 3)+
                      theme(text = element_text(size=18),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                            axis.title = element_text(size = 15),
                            legend.title = element_blank())+
                      geom_vline(xintercept = 0, linetype = "dashed")+
                      geom_hline(yintercept = 0, linetype = "dashed")+
                      geom_text(aes(label=Population),hjust=-0.1, vjust=1)+
                      ylim( -6.2, 1)+
                      xlim(-0.85, 0.3)+
                      ylab(bquote(Delta~Establishment~Rate))+
                      xlab(bquote(Delta~PC1))

delta_PC1_BRTE <- ggplot(delta.PCA.survival%>% filter(Competition == "BRTE"), aes(x = delta_PC1, y = delta_POSE))+
                      geom_point(colour = "#777777", size = 3)+
                      theme(text = element_text(size=18),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                            axis.title = element_text(size = 15),
                            legend.title = element_blank())+
                      geom_vline(xintercept = 0, linetype = "dashed")+
                      geom_hline(yintercept = 0, linetype = "dashed")+
                      geom_text(aes(label=Population),hjust=-0.1, vjust=1)+
                      ylim( -6.2, 1)+
                      xlim(-0.85, 0.3)+
                      ylab(bquote(Delta~Establishment~Rate))+
                      xlab(bquote(Delta~PC1))

delta_PC2_None <- ggplot(delta.PCA.survival%>% filter(Competition == "None"), aes(x = delta_PC2, y = delta_POSE))+
                      geom_point(colour = "#111111", size = 3)+
                      theme(text = element_text(size=18),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                            axis.title = element_text(size = 15),
                            legend.title = element_blank())+
                      geom_vline(xintercept = 0, linetype = "dashed")+
                      geom_hline(yintercept = 0, linetype = "dashed")+
                      geom_text(aes(label=Population),hjust=-0.1, vjust=1)+
                      ylim( -6, 1)+
                      xlim(-0.1, 0.62)+
                      ylab(bquote(Delta~Establishment~Rate))+
                      xlab(bquote(Delta~PC2))

delta_PC2_BRTE <- ggplot(delta.PCA.survival%>% filter(Competition == "BRTE"), aes(x = delta_PC2, y = delta_POSE))+
                      geom_point(colour = "#777777", size = 3)+
                      theme(text = element_text(size=18),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                            axis.title = element_text(size = 15),
                            legend.title = element_blank())+
                      geom_vline(xintercept = 0, linetype = "dashed")+
                      geom_hline(yintercept = 0, linetype = "dashed")+
                      geom_text(aes(label=Population),hjust=0.9, vjust=1.2)+
                      ylim( -6, 1)+
                      xlim(-0.1, 0.62)+
                      ylab(bquote(Delta~Establishment~Rate))+
                      xlab(bquote(Delta~PC2))
ggarrange(delta_PC1_None, delta_PC2_None, distance_None, delta_PC1_BRTE, delta_PC2_BRTE, distance_BRTE,
          ncol = 3, nrow = 2, labels = c("None","None", "None", "BRTE", "BRTE", "BRTE"), 
          label.x = .1, label.y = 0.99)

#----------------------------------------------------------#
# Is trait plasticity correlated with drought tolerance?
# Delta trait (Dry-Wet) vs delta establishment rate (Dry-Wet) for each trait
mean_survival_trait <- survival_trait %>%
  dplyr::select(-TotalBiomass)%>%
  pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value")%>%
  filter(trait%in%c("Emergence", "Height", "SLA", "LDMC", "RMR",  "Tips", "PropF", "SRL", "Length")) %>%
  group_by(Population, Water, Competition, Treatment, trait) %>%
  summarise(mean_trait = mean(value),
            mean_POSE = mean(POSE_survival_stem_count/25))

delta_survival_trait <- mean_survival_trait %>%
  group_by(Population, Competition, trait) %>%
  summarise(delta_trait = abs((mean_trait[Water == "Dry"]-mean_trait[Water == "Wet"])/mean_trait[Water == "Wet"]),
            delta_est = (mean_POSE[Water == "Dry"]-mean_POSE[Water == "Wet"]))

summary(lm(delta_est ~ delta_trait, data = delta_survival_trait%>%filter(trait == "Emergence"))) #*
summary(lm(delta_est ~ delta_trait, data = delta_survival_trait%>%filter(trait == "Height"))) #NA
summary(lm(delta_est ~ delta_trait, data = delta_survival_trait%>%filter(trait == "LDMC"))) #NA
summary(lm(delta_est ~ delta_trait, data = delta_survival_trait%>%filter(trait == "Length"))) #*
summary(lm(delta_est ~ delta_trait, data = delta_survival_trait%>%filter(trait == "PropF"))) #NA
summary(lm(delta_est ~ delta_trait, data = delta_survival_trait%>%filter(trait == "RMR"))) #NA
summary(lm(delta_est ~ delta_trait, data = delta_survival_trait%>%filter(trait == "SLA"))) #NA
summary(lm(delta_est ~ delta_trait, data = delta_survival_trait%>%filter(trait == "SRL"))) #NA
summary(lm(delta_est ~ delta_trait, data = delta_survival_trait%>%filter(trait == "Tips"))) #NA

ann_text_plastic <- data.frame(delta_trait = c(1.4, 0.6),
                               delta_est = 0.02,
                               label = c("R2 = 0.36, p = 0.03", "R2 = 0.37, p = 0.03"),
                       trait =  c("Emergence","Length"))
library(plyr)
delta_survival_trait$Population <- revalue(delta_survival_trait$Population, c("Butte Valley"="B", "Steens"="S", "EOARC"="E", "Water Canyon"="W", "Reno"="R"))

ggplot(delta_survival_trait, aes(x = delta_trait, y = delta_est))+
  geom_point(aes(col = Competition), size = 3)+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15))+
  facet_wrap(~trait, scale = "free", ncol = 3, labeller=trait_labeller)+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_text(aes(col = Competition, label=Population),hjust="inward", vjust="inward",show.legend = FALSE)+
  #ylim( -1, 0.1)+
  #xlim(0, 1.1)+
  ylab(bquote(Drought~Tolerance))+
  xlab(bquote(Trait~Plasticity))+
  scale_color_manual(values=c("#FF9507","#CE026E"))+
  geom_smooth(data = delta_survival_trait%>%filter(trait %in% c("Emergence", "Length")), col = "black",method = "lm",  size=0.9, se = FALSE)+
  geom_text( data = ann_text_plastic, mapping = aes(x = delta_trait, y = delta_est, label = label))
# # Calculate PCA score distance as a metric of plasticity
# plasticity <- pca_trait_scores_lab %>%
#   dplyr::select(Population, Replicate, Competition, Water, PC1, PC2) %>%
#   pivot_wider(names_from = Water, values_from = c(PC1, PC2)) %>%
#   group_by(Population, Replicate, Competition) %>%
#   na.omit()%>%
#   summarise(d_trait = sqrt((PC1_Dry-PC1_Wet)^2 + (PC2_Dry- PC2_Wet)^2))
# 
# # Calculate drought tolerance: relative change in survival rate and biomass 
# survival_delta <- growth%>%
#   dplyr::select(Population, Replicate, Water, Competition, POSE_survival_stem_count) %>%
#   mutate(survivalrate = POSE_survival_stem_count/25) %>%
#   dplyr::select(Population, Replicate, Water, Competition, survivalrate) %>%
#   pivot_wider(names_from = Water, values_from = survivalrate) %>%
#   group_by(Population, Replicate, Competition) %>%
#   summarise(survival_delta = (Dry-Wet)/Wet)
# biomass_delta <- biomass_dat %>%
#   filter(Life_stage == "seedling") %>%
#   dplyr::select(Population, Replicate, Water, Competition, TotalBiomass) %>%
#   pivot_wider(names_from = Water, values_from = TotalBiomass) %>%
#   group_by(Population, Replicate, Competition) %>%
#   na.omit()%>%
#   summarise(biomass_delta = (Dry-Wet)/Wet)
# 
# # Combine dataframes to plot it
# plastic_combined <- left_join(survival_delta, biomass_delta) %>%
#   left_join(., plasticity) 
# f1 <- ggplot(plastic_combined, aes(x = d_trait, y = survival_delta, col = Competition)) +
#   geom_point(aes(col = Competition))+
#   theme(text = element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#         axis.title = element_text(size = 12))+
#   ylab(bquote(Relative~Change~Establishment~Rate)) +
#   xlab("Trait plasticity")+
#   geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5)
# f2 <- ggplot(plastic_combined, aes(x = d_trait, y = biomass_delta, col = Competition)) +
#   geom_point(aes(col = Competition))+
#   theme(text = element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#         axis.title = element_text(size = 12))+
#   ylab(bquote(Relative~Change~Total~Biomass)) +
#   xlab("Trait plasticity")+
#   geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5)
# 
# 
# # Graph them together
# ggarrange(f1, f2, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
#           font.label = list(size = 15), common.legend = TRUE, legend = "right", label.x = -.02, label.y = .99)
# 
# # Alternative:
# # Calculate trait plasticity or deviation from None_Wet as control 
# plasticity_v2 <- pca_trait_scores_lab %>%
#   dplyr::select(Population, Replicate, Treatment, PC1, PC2) %>%
#   pivot_wider(names_from = Treatment, values_from = c(PC1, PC2)) %>%
#   group_by(Population, Replicate) %>%
#   na.omit()%>%
#   summarise(None_Dry = sqrt((PC1_None_Dry-PC1_None_Wet)^2 + (PC2_None_Dry- PC2_None_Wet)^2),
#             BRTE_Wet = sqrt((PC1_BRTE_Wet-PC1_None_Wet)^2 + (PC2_BRTE_Wet- PC2_None_Wet)^2),
#             BRTE_Dry = sqrt((PC1_BRTE_Dry-PC1_None_Wet)^2 + (PC2_BRTE_Dry- PC2_None_Wet)^2)) %>%
#   pivot_longer(cols = None_Dry:BRTE_Dry, names_to = "Treatment", values_to = "d_trait")
# 
# ggplot(plasticity_v2, aes(x = Population, y = d_trait, fill = Treatment)) +
#   geom_boxplot()+
#   theme_classic()+
#   scale_fill_manual(values=c("#999999", "#6A0DAD", "#E69F00"))+
#   ylab(bquote(Trait~Plasticity)) 
# 
# 
# # Calculate stress tolerance: relative change in survival rate and biomass 
# growth$Treatment <- apply(growth[ ,3:4 ] , 1 , paste , collapse = "_" )
# biomass_dat$Treatment <- apply(biomass_dat[ ,3:4 ] , 1 , paste , collapse = "_" )
# survival_delta_v2 <- growth%>%
#   dplyr::select(Population, Replicate, Treatment, POSE_survival_stem_count) %>%
#   mutate(survivalrate = POSE_survival_stem_count/25) %>%
#   dplyr::select(Population, Replicate, Treatment, survivalrate) %>%
#   pivot_wider(names_from = Treatment, values_from = survivalrate) %>%
#   group_by(Population, Replicate) %>%
#   summarise(None_Dry = (None_Dry-None_Wet)/None_Wet,
#             BRTE_Wet = (BRTE_Wet-None_Wet)/None_Wet,
#             BRTE_Dry = (BRTE_Dry-None_Wet)/None_Wet) %>%
#   pivot_longer(cols = None_Dry:BRTE_Dry, names_to = "Treatment", values_to = "survival_delta")
# biomass_delta_v2 <- biomass_dat %>%
#   filter(Life_stage == "seedling") %>%
#   dplyr::select(Population, Replicate, Treatment, TotalBiomass) %>%
#   pivot_wider(names_from = Treatment, values_from = TotalBiomass) %>%
#   group_by(Population, Replicate) %>%
#   na.omit()%>%
#   summarise(None_Dry = (None_Dry-None_Wet)/None_Wet,
#             BRTE_Wet = (BRTE_Wet-None_Wet)/None_Wet,
#             BRTE_Dry = (BRTE_Dry-None_Wet)/None_Wet) %>%
#   pivot_longer(cols = None_Dry:BRTE_Dry, names_to = "Treatment", values_to = "biomass_delta")
# 
# # Combine dataframes to plot it
# plastic_combined_v2 <- left_join(survival_delta_v2, biomass_delta_v2) %>%
#   left_join(., plasticity_v2) 
# f1_V2 <- ggplot(plastic_combined_v2, aes(x = d_trait, y = survival_delta)) +
#   geom_point(aes(col = Treatment))+facet_grid(cols = vars(Population))+
#   theme(text = element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#         axis.title = element_text(size = 12))+
#   ylab(bquote(Relative~Change~Establishment~Rate)) +
#   xlab("Trait plasticity")+
#   geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5, color = "black")+
#   scale_color_manual(values=c("#999999", "#6A0DAD", "#E69F00"))
# f2_v2 <- ggplot(plastic_combined_v2, aes(x = d_trait, y = biomass_delta)) +
#   geom_point(aes(col = Treatment))+facet_grid(cols = vars(Population))+
#   theme(text = element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#         axis.title = element_text(size = 12))+
#   ylab(bquote(Relative~Change~Total~Biomass)) +
#   xlab("Trait plasticity")+
#   geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5, color = "black")+
#   scale_color_manual(values=c("#999999", "#6A0DAD", "#E69F00"))
# 
# # Graph them together
# ggarrange(f1_V2, f2_v2, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
#           font.label = list(size = 15), common.legend = TRUE, legend = "right", label.x = -.02, label.y = .99)
# 
# # Plasticity of PropF and LDMC and relative change in establishment rate and biomass
# trait_master$Treatment <- apply(trait_master[ ,3:4 ] , 1 , paste , collapse = "_" )
# PropF_delta <- trait_master %>%
#   dplyr::select(Population, Replicate, Treatment, PropF) %>%
#   pivot_wider(names_from = Treatment, values_from = PropF) %>%
#   group_by(Population, Replicate) %>%
#   summarise(None_Dry = (None_Dry-None_Wet)/None_Wet,
#             BRTE_Wet = (BRTE_Wet-None_Wet)/None_Wet,
#             BRTE_Dry = (BRTE_Dry-None_Wet)/None_Wet) %>%
#   pivot_longer(cols = None_Dry:BRTE_Dry, names_to = "Treatment", values_to = "PropF_delta")
# 
# LDMC_delta <- trait_master %>%
#   dplyr::select(Population, Replicate, Treatment, LDMC) %>%
#   pivot_wider(names_from = Treatment, values_from = LDMC) %>%
#   group_by(Population, Replicate) %>%
#   summarise(None_Dry = (None_Dry-None_Wet)/None_Wet,
#             BRTE_Wet = (BRTE_Wet-None_Wet)/None_Wet,
#             BRTE_Dry = (BRTE_Dry-None_Wet)/None_Wet) %>%
#   pivot_longer(cols = None_Dry:BRTE_Dry, names_to = "Treatment", values_to = "LDMC_delta")
# 
# trait_delta_combined <- left_join(survival_delta_v2, biomass_delta_v2) %>%
#   left_join(., PropF_delta)%>%
#   left_join(., LDMC_delta)
# trait_delta_combined$Population <- ordered(as.factor(trait_delta_combined$Population), levels = c("Butte Valley","Steens","EOARC",
#                                                                                   "Water Canyon",  "Reno"))
# f1_V3 <- ggplot(trait_delta_combined, aes(x = PropF_delta, y = survival_delta)) +
#   geom_jitter(aes(col = Treatment))+
#   facet_grid(cols = vars(Population))+
#   theme(text = element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#         axis.title = element_text(size = 12))+
#   ylab(bquote(Relative~Change~Establishment~Rate)) +
#   xlab("PropF plasticity")+
#   geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5, color = "black")+
#   scale_color_manual(values=c("#999999", "#6A0DAD", "#E69F00"))
# f2_v3 <- ggplot(trait_delta_combined, aes(x = PropF_delta, y = biomass_delta)) +
#   geom_point(aes(col = Treatment))+
#   facet_grid(cols = vars(Population))+
#   theme(text = element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#         axis.title = element_text(size = 12))+
#   ylab(bquote(Relative~Change~Total~Biomass)) +
#   xlab("PropF plasticity")+
#   geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5, color = "black")+
#   scale_color_manual(values=c("#999999", "#6A0DAD", "#E69F00"))
# f3_v3 <- ggplot(trait_delta_combined, aes(x = LDMC_delta, y = survival_delta)) +
#   geom_jitter(aes(col = Treatment))+
#   facet_grid(cols = vars(Population))+
#   theme(text = element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#         axis.title = element_text(size = 12))+
#   ylab(bquote(Relative~Change~Establishment~Rate)) +
#   xlab("LDMC plasticity")+
#   geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5, color = "black")+
#   scale_color_manual(values=c("#999999", "#6A0DAD", "#E69F00"))
# f4_v3 <- ggplot(trait_delta_combined, aes(x = LDMC_delta, y = biomass_delta)) +
#   geom_point(aes(col = Treatment))+
#   facet_grid(cols = vars(Population))+
#   theme(text = element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#         axis.title = element_text(size = 12))+
#   ylab(bquote(Relative~Change~Total~Biomass)) +
#   xlab("LDMC plasticity")+
#   geom_smooth(method = "lm", formula = y~poly(x,2), size=0.5, color = "black")+
#   scale_color_manual(values=c("#999999", "#6A0DAD", "#E69F00"))
# # Graph them together
# ggarrange(f1_V3, f2_v3, f3_v3, f4_v3, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"),
#           font.label = list(size = 15), common.legend = TRUE, legend = "right", label.x = -.02, label.y = .99)

#--Seed mass summary----#
summarizeseeds <- seeds %>% 
  group_by(Population,Species) %>%
  filter(!is.na(Weight_g)) %>%
  filter(Population %in% c("Butte Valley", "Steens", "EOARC", "Water Canyon")) %>%
  summarise(meanseed=mean(Weight_g),seseed=se(Weight_g))
summarizeseeds$Population <- ordered(as.factor(summarizeseeds$Population), levels = c("Butte Valley","Steens","EOARC",
                                                                                                  "Water Canyon"))
summary(aov(Weight_g~ Population, seeds))
f_seeds <- ggplot(summarizeseeds, aes(x = Population, y = meanseed/10))+
                geom_point()+
                theme(text = element_text(size=15),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position = "none",
                      legend.title = element_blank(),
                      axis.title = element_text(size = 15))+
                geom_errorbar(aes(ymin = meanseed/10-seseed/10, ymax = meanseed/10+seseed/10), 
                              width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                ylab("Seed mass (g)")
# Did seed mass relate to any seedling traits??
seedsmass_to_traits <- trait_master %>%
  inner_join(., summarizeseeds)
pairs(~meanseed + Emergence+ PropF+Length, seedsmass_to_traits)
# Did seed mass relate to establishment rate?
seedsmass_to_growth <- growth %>%
  inner_join(., summarizeseeds)
pairs(~meanseed + BRTE_stem_count + POSE_survival_stem_count, seedsmass_to_growth)
f_seed_establishment <- ggplot(seedsmass_to_growth, aes(x = meanseed/10, y = POSE_survival_stem_count/25))+
                            geom_jitter()+
                            theme(text = element_text(size=15),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.background = element_blank(),
                                  axis.line = element_line(colour = "black"),
                                  legend.position = "none",
                                  legend.title = element_blank(),
                                  axis.title = element_text(size = 15))+
                                  #axis.title.y = element_blank())
                                  #axis.title.x = element_blank(),
                                  #axis.text.x = element_text(angle = 90))
                            geom_smooth(method = lm)+
                            ylab("Establishment Rate")+
                            xlab("Seed mass (g)")
summary(lm(POSE_survival_stem_count ~ meanseed, seedsmass_to_growth))
# Did seed mass relate to biomass?
f_seed_biomass <- ggplot(seedsmass_to_traits, aes(x = meanseed/10, y = TotalBiomass))+
                      geom_jitter()+
                      geom_smooth(method = lm)+
                      theme(text = element_text(size=15),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            legend.position = "none",
                            legend.title = element_blank(),
                            axis.title = element_text(size = 15))+
                      ylab("Total biomass (g)")+
                      xlab("Seed mass (g)")
summary(lm(TotalBiomass ~ meanseed, seedsmass_to_traits))
ggarrange(f_seeds, f_seed_establishment, f_seed_biomass, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"),
                   font.label = list(size = 12))
