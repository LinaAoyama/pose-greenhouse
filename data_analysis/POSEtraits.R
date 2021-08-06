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

trait_master <- inner_join(potID, root_traits)%>%
  inner_join(.,biomass_traits) %>%
  inner_join(., height_trait) %>%
  inner_join(., leaf_traits)

trait_matrix <- as.matrix(trait_master[,7:ncol(trait_master)]) #extract only the trait value columns

#----------------------------------------------------------#
# What are the trait characteristics of each population?
# Create NMDS of traits
nmds_trait = metaMDS(trait_matrix, distance = "bray") #generate an NMDS plot
plot(nmds_trait) #visualize the plot
nmds_trait_scores = as.data.frame(scores(nmds_trait)) #extract nmds scores
nmds_trait_scores_lab = as.data.frame(cbind(trait_master[,1:6],nmds_trait_scores)) #add plot info back
nmds_trait_scores_lab$Population <- ordered(as.factor(nmds_trait_scores_lab$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                           "Water Canyon",  "Reno"))
en = envfit(nmds_trait, trait_matrix, permutations = 999, na.rm = TRUE)
en_coord_cont = as.data.frame(as.data.frame(scores(en, "vectors"))* ordiArrowMul(en))
NMDS_trait <- ggplot(nmds_trait_scores_lab, aes(x = NMDS1, y = NMDS2)) + #plot nmds in ggplot2
                  geom_point(size = 4, aes( colour = Population), alpha = 0.5)+
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
                        axis.title = element_text(size = 15))+
                  geom_segment(data = en_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                            alpha = 0.5, size = 1, colour = "grey30") +
                  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                            fontface = "bold", label = row.names(en_coord_cont)) +
                  xlim(-0.8, 0.4)

# Plot spread of trait values by population
trait_long <- trait_master %>%
  pivot_longer(cols = Length_cm:LDMC, names_to = "trait", values_to = "value")
trait_long$Population <- ordered(as.factor(trait_long$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                                                          "Water Canyon",  "Reno"))
trait_long$Treatment <- apply(trait_long[ ,3:4 ] , 1 , paste , collapse = "-" )
ggplot(trait_long, aes(x = Population, y = value, col = Competition)) +
  geom_violin()+
  theme_bw()+
  #geom_point()+
  facet_wrap(~trait, scales = "free", ncol = 2)
ggplot(trait_long, aes(x = Population, y = value, col = Water)) +
  geom_violin()+
  #geom_point()+
  facet_wrap(~trait, scales = "free", ncol = 2)
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
#----------------------------------------------------------#
# How did the trait shift by BRTE competition and water treatment?
# Create separate NMDS for aboveground and belowground traits

# Aboveground traits
AG_trait_matrix <- as.matrix(trait_master[,15:ncol(trait_master)])
AG_nmds_trait = metaMDS(AG_trait_matrix, distance = "bray")
AG_nmds_trait_scores = as.data.frame(scores(AG_nmds_trait))
AG_nmds_trait_scores_lab = as.data.frame(cbind(trait_master[,1:6],AG_nmds_trait_scores)) 
AG_nmds_trait_scores_lab$Population <- ordered(as.factor(AG_nmds_trait_scores_lab$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                                                    "Water Canyon",  "Reno"))
AG_nmds_trait_scores_lab$Treatment <- apply(AG_nmds_trait_scores_lab[ ,3:4 ] , 1 , paste , collapse = "-" )
AG_en = envfit(AG_nmds_trait, AG_trait_matrix, permutations = 999, na.rm = TRUE)
AG_en_coord_cont = as.data.frame(as.data.frame(scores(AG_en, "vectors"))* ordiArrowMul(AG_en))
AG_NMDS_trait <- ggplot(AG_nmds_trait_scores_lab, aes(x = NMDS1, y = NMDS2)) + #plot nmds in ggplot2
                  geom_point(size = 4, aes( shape = Competition, colour = Treatment), alpha = 0.5)+
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
                        axis.title = element_text(size = 15))+
                  geom_segment(data = AG_en_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                               alpha = 0.5, size = 1, colour = "grey30") +
                  geom_text(data = AG_en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                            fontface = "bold", label = row.names(AG_en_coord_cont)) +
                  xlim(-0.8, 0.5)

# Belowground traits
BG_trait_matrix <- as.matrix(trait_master[,7:14])
BG_nmds_trait = metaMDS(BG_trait_matrix, distance = "bray")
BG_nmds_trait_scores = as.data.frame(scores(BG_nmds_trait))
BG_nmds_trait_scores_lab = as.data.frame(cbind(trait_master[,1:6],BG_nmds_trait_scores)) 
BG_nmds_trait_scores_lab$Population <- ordered(as.factor(BG_nmds_trait_scores_lab$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                                                          "Water Canyon",  "Reno"))
BG_nmds_trait_scores_lab$Treatment <- apply(BG_nmds_trait_scores_lab[ ,3:4 ] , 1 , paste , collapse = "-" )
BG_en = envfit(BG_nmds_trait, BG_trait_matrix, permutations = 999, na.rm = TRUE)
BG_en_coord_cont = as.data.frame(as.data.frame(scores(BG_en, "vectors"))* ordiArrowMul(BG_en))
BG_NMDS_trait <- ggplot(BG_nmds_trait_scores_lab, aes(x = NMDS1, y = NMDS2)) + #plot nmds in ggplot2
                geom_point(size = 4, aes( shape = Competition, colour = Treatment), alpha = 0.5)+
                theme(text = element_text(size=15),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
                      axis.title = element_text(size = 15))+
                geom_segment(data = BG_en_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                             alpha = 0.5, size = 1, colour = "grey30") +
                geom_text(data = BG_en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                          fontface = "bold", label = row.names(BG_en_coord_cont), 
                          position=position_jitter(width=-0.5,height=0)) +
                xlim(-0.8, 0.7)

# Graph them together
ggarrange(AG_NMDS_trait, BG_NMDS_trait, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
          font.label = list(size = 15), common.legend = TRUE, legend = "right")
