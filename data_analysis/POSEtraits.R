# Preliminary data exploration
# Seedling traits of Sandberg bluegrass (POSE) in the greenhouse

# Set pathway first

# Data
source("data_compiling/compile_demography.R")

# Packages
library(tidyverse) #data wrangling
library(ggplot2) #plot
library(ggpubr) #combine plots
library(vegan) #nmds
library(corrplot) #correlation matrix

# Prep data
root_traits <- full_join(rootbiomass, root) %>%
  group_by(PotID) %>%
  mutate(SRL = Length_cm/Root_Weight_g,
         Coarse_cm = LenTotHistoClasses-D0L05,
         Fine_cm = D0L05) %>% #calculate specific root length, coarse root length, and fine root length
  select(PotID, Length_cm, SRL, SurfArea_cm2, AvgDiam_mm, Tips, Forks, Fine_cm, Coarse_cm) 
  
biomass_traits <- biomass %>%
  filter(Species == "POSE") %>%
  inner_join(., rootbiomass) %>%
  select(-Species) %>%
  mutate(TotalBiomass = Dry_Biomass_Weight_g+Root_Weight_g,
         RMR = Root_Weight_g/TotalBiomass) %>% #calculate total biomass weight and root to mass ratio
  select(PotID, TotalBiomass, RMR)

height_trait<- height %>% 
  group_by(PotID) %>% 
  summarise(Height = mean(Height_cm)) #calculate average seedling height

leaf_traits <- leaftraits %>%
  drop_na() %>%
  group_by(PotID) %>%
  summarise(SLA = FreshLeafArea_cm2/DryLeafWeight_g/NumberLeaves, 
         LDMC = DryLeafWeight_g/FreshLeafArea_cm2/NumberLeaves) #calculate specific leaf area and leaf dry matter content

germination_dates <- germination %>%
  select(PotID, Days_emergence) %>%
  drop_na()

trait_master <- inner_join(potID, root_traits)%>%
  inner_join(.,biomass_traits) %>%
  inner_join(., height_trait) %>%
  inner_join(., leaf_traits) %>%
  inner_join(., germination_dates) %>%
  filter(Population != "Gund")

trait_matrix <- as.matrix(trait_master[,7:ncol(trait_master)]) #extract only the trait value columns

# Function for standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
} 

#----------------------------------------------------------#
# What are the trait characteristics of each population?
# Create PCA of traits
pca_trait = rda(trait_matrix, scale = TRUE) #run PCA on all traits
biplot(pca_trait, display = c("sites", "species"), type = c("text", "points")) #plot biplot
pca_trait_scores <- as.data.frame(scores(pca_trait, choices=c(1,2), display=c("sites"))) #extract pca1 and pca2 scores
pca_trait_scores_lab = as.data.frame(cbind(trait_master[,1:6],pca_trait_scores)) #add plot info back
pca_trait_scores_lab$Population <- ordered(as.factor(pca_trait_scores_lab$Population), levels = c("Butte Valley","Steens","EOARC", 
                                                                                                    "Water Canyon",  "Reno"))
envout<-as.data.frame(scores(pca_trait, choices=c(1,2), display=c("species")))
ggplot(pca_trait_scores_lab, aes(x = PC1, y = PC2))+
  geom_point(size = 4, aes(colour = Population), alpha = 0.5)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 15))+
  geom_segment(data = envout, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text(data = envout, aes(x = PC1, y = PC2), colour = "grey30", 
            fontface = "bold", label = row.names(envout))+
  xlim(-2, 2.5)

# Plot spread of trait values by population
trait_long <- trait_master %>%
  pivot_longer(cols = Length_cm:Days_emergence, names_to = "trait", values_to = "value")
trait_long$Population <- ordered(as.factor(trait_long$Population), levels = c("Butte Valley","Steens","EOARC", 
                                                                              "Water Canyon",  "Reno"))
trait_long$Treatment <- apply(trait_long[ ,3:4 ] , 1 , paste , collapse = "-" )
mean_trait_long <- trait_long%>%
  group_by(Population, Treatment, trait) %>%
  summarise(mean = mean(value),
            se = se(value))
f_AG_trait <-ggplot(mean_trait_long%>%filter(Treatment == "None-Wet")%>%filter(trait%in%c("Days_emergence", "Height", "SLA", "LDMC")), aes(x = Population, y = mean)) +
                geom_point()+
                geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                facet_wrap(~trait, scales = "free", ncol = 2)+
                theme_bw()+
                ylab("")
f_BG_trait <-ggplot(mean_trait_long%>%filter(Treatment == "None-Wet")%>%filter(!(trait%in%c("Days_emergence", "Height", "SLA", "LDMC"))), aes(x = Population, y = mean)) +
                geom_point()+
                geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                facet_wrap(~trait, scales = "free", ncol = 2)+
                theme_bw()+
                ylab("")

# Graph them together
ggarrange(f_AG_trait, f_BG_trait, ncol = 1, nrow = 2, labels = c("(a)", "(b)"),
          font.label = list(size = 15), common.legend = TRUE, legend = "right", heights = c(1, 2.5))
TukeyHSD(aov(value~Population, data = trait_long%>%
               filter(Treatment == "None-Wet")%>%
               filter(trait == "SLA")))

#----------------------------------------------------------#
# How did the trait shift by BRTE competition and water treatment?
# Create separate PCAs for aboveground and belowground traits

# Aboveground traits
AG_trait_matrix <- as.matrix(trait_master[,17:ncol(trait_master)])
AG_pca_trait = rda(AG_trait_matrix, scale = TRUE)
AG_pca_trait_scores = as.data.frame(scores(AG_pca_trait, choices=c(1,2), display=c("sites")))
AG_pca_trait_scores_lab = as.data.frame(cbind(trait_master[,1:6],AG_pca_trait_scores)) 
AG_pca_trait_scores_lab$Population <- ordered(as.factor(AG_nmds_trait_scores_lab$Population), levels = c("Butte Valley","Steens","EOARC",
                                                                                                          "Water Canyon",  "Reno"))
AG_pca_trait_scores_lab$Treatment <- apply(AG_pca_trait_scores_lab[ ,3:4 ] , 1 , paste , collapse = "-" )
summary(AG_pca_trait)
AG_en = as.data.frame(scores(AG_pca_trait, choices=c(1,2), display=c("species")))
f_AG_pca <- ggplot(AG_pca_trait_scores_lab, aes(x = PC1, y = PC2)) + #plot pca in ggplot2
  geom_point(size = 4, aes( shape = Population, colour = Treatment), alpha = 0.5)+
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 15))+
  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
  geom_segment(data = AG_en, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text(data = AG_en, aes(x = PC1, y = PC2), colour = "grey30", 
            fontface = "bold", label = row.names(AG_en))+
  xlab("PC1 (37.2%)")+
  ylab("PC2 (29.7%)")+
  xlim(-2.3,2.3)

# Belowground traits
BG_trait_matrix <- as.matrix(trait_master[,7:16])
BG_pca_trait = rda(BG_trait_matrix, scale = TRUE)
BG_pca_trait_scores = as.data.frame(scores(BG_pca_trait, choices = c(1,2), display = c("sites")))
BG_pca_trait_scores_lab = as.data.frame(cbind(trait_master[,1:6],BG_pca_trait_scores)) 
BG_pca_trait_scores_lab$Population <- ordered(as.factor(BG_pca_trait_scores_lab$Population), levels = c("Butte Valley","Steens","EOARC", 
                                                                                                          "Water Canyon",  "Reno"))
BG_pca_trait_scores_lab$Treatment <- apply(BG_pca_trait_scores_lab[ ,3:4 ] , 1 , paste , collapse = "-" )
summary(BG_pca_trait)
BG_en = as.data.frame(scores(BG_pca_trait, choices=c(1,2), display=c("species")))
f_BG_pca <- ggplot(BG_pca_trait_scores_lab, aes(x = PC1, y = PC2)) + #plot pca in ggplot2
  geom_point(size = 4, aes( shape = Population, colour = Treatment), alpha = 0.5)+
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        axis.title = element_text(size = 15))+
  scale_shape_manual(values=c(8, 5, 15, 17, 19))+
  geom_segment(data = BG_en, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               alpha = 0.5, size = 1, colour = "grey30") +
  geom_text(data = BG_en, aes(x = PC1, y = PC2), colour = "grey30", 
            fontface = "bold", label = row.names(BG_en))+
  xlab("PC1 (70.2%)")+
  ylab("PC2 (14.6%)")+
  ylim(-1.5, 1.2)+
  xlim(-1.5,2.3)

# Graph them together
ggarrange(f_AG_pca, f_BG_pca, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
          font.label = list(size = 15), common.legend = TRUE, legend = "right")


#----------------------------------------------------------#
# Which trait is correlated with higher POSE survival rate?
survival_trait <- inner_join(trait_master, stemcount) %>%
  select(-POSE_emergence_stem_count) %>%
  select(-BRTE_stem_count) %>%
  drop_na()
  
# This function calculates p-value of the correlation
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(survival_trait[,7:ncol(survival_trait)])# matrix of the p-values

corrplot(cor(survival_trait[,7:ncol(survival_trait)]), type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.01, insig = "blank")

ggplot(survival_trait%>%pivot_longer(cols = Length_cm:LDMC, names_to = "trait", values_to = "value"), aes(x = value, y = POSE_survival_stem_count))+
  geom_point()+
  facet_wrap(~trait, scale = "free")
