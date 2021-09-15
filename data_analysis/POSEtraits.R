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
colnames(trait_master) <- c('PlotID', 'Life_stage', 'Competition', 'Water', 'Population', 'Replicate',
                            'Length', 'SRL', 'SurfArea', 'AvgDiam', 'Tips', 'Forks', 'Fine', 'Coarse', 
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
# Take out any traits that are collinear
# Standardize the data
pairs(~RMR + Tips + Length + Fine + Coarse + SRL + SurfArea + AvgDiam + Forks + TotalBiomass, trait_matrix_raw)
trait_matrix <- as.data.frame(decostand(trait_matrix_raw, "standardize")) %>%
  dplyr::select( -Fine, -Forks, -SRL, -Coarse, -SurfArea, -TotalBiomass)


#----------------------------------------------------------#
# How did the trait shift by BRTE competition and water treatment?
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
  geom_point(size = 4, aes(colour = Treatment, shape = Population), alpha = 0.5)+
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
  xlim(-2, 2.3)

#----------------------------------------------------------#
# How did each trait respond to the treatment?
# Calculate mean and standard error of each trait by population and treatment
trait_long <- trait_master %>%
  pivot_longer(cols = Length:Emergence, names_to = "trait", values_to = "value")
trait_long$Population <- ordered(as.factor(trait_long$Population), levels = c("Butte Valley","Steens","EOARC", 
                                                                              "Water Canyon",  "Reno"))
trait_long$Treatment <- apply(trait_long[ ,3:4 ] , 1 , paste , collapse = "-" )
mean_trait_long <- trait_long%>%
  group_by(Population, Treatment, trait) %>%
  summarise(mean = mean(value),
            se = se(value))
f_AG_trait <-ggplot(mean_trait_long%>%filter(trait%in%c("Emergence", "Height", "SLA", "LDMC")), aes(x = Treatment, y = mean, col = Population)) +
                geom_point(position = position_dodge(width = 0.5))+
                geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                facet_wrap(~trait, scales = "free", ncol = 2)+
                theme_bw()+
                ylab("")
f_BG_trait <-ggplot(mean_trait_long%>%filter((trait%in%c("AvgDiam", "RMR", "Length", "Tips"))), aes(x = Treatment, y = mean, col = Population)) +
                geom_point(position = position_dodge(width = 0.5))+
                geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                facet_wrap(~trait, scales = "free", ncol = 2)+
                theme_bw()+
                ylab("")

# Graph them together
ggarrange(f_AG_trait, f_BG_trait, ncol = 1, nrow = 2, labels = c("(a)", "(b)"),
          font.label = list(size = 15), common.legend = TRUE, legend = "right", heights = c(1, 1))



#----------------------------------------------------------#
# Is trait plasticity correlated with BRTE resistance?

# Calculate PCA score distance as a metric of plasticity
plasticity <- pca_trait_scores_lab %>%
  dplyr::select(Population, Replicate, Competition, Water, PC1, PC2) %>%
  pivot_wider(names_from = Competition, values_from = c(PC1, PC2)) %>%
  group_by(Population, Replicate, Water) %>%
  na.omit()%>%
  summarise(d_trait = sqrt((PC1_BRTE-PC1_None)^2 + (PC2_BRTE- PC2_None)^2))

# Calculate BRTE resistance: survival rate BRTE - None and biomass BRTE - None
survival_delta <- growth%>%
  dplyr::select(Population, Replicate, Water, Competition, POSE_survival_stem_count) %>%
  mutate(survivalrate = POSE_survival_stem_count/25) %>%
  dplyr::select(Population, Replicate, Water, Competition, survivalrate) %>%
  pivot_wider(names_from = Competition, values_from = survivalrate) %>%
  group_by(Population, Replicate, Water) %>%
  summarise(survival_delta = BRTE-None)
biomass_delta <- biomass_dat %>%
  filter(Life_stage == "seedling") %>%
  dplyr::select(Population, Replicate, Water, Competition, POSE) %>%
  pivot_wider(names_from = Competition, values_from = POSE) %>%
  group_by(Population, Replicate, Water) %>%
  na.omit()%>%
  summarise(biomass_delta = BRTE-None)

# Combine dataframes to plot it
plastic_combined <- left_join(survival_delta, biomass_delta) %>%
  left_join(., plasticity) 
f1 <- ggplot(plastic_combined, aes(x = d_trait, y = survival_delta)) +
        geom_point(aes(col = Water))+
        theme(text = element_text(size=12),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
              axis.title = element_text(size = 12))+
        ylab(bquote(italic(P.~secunda)~Survival~Rate~(BRTE-None))) +
        xlab("Trait platicity")+
        geom_smooth(method = "lm", colour="black", size=0.5)
f2 <- ggplot(plastic_combined, aes(x = d_trait, y = biomass_delta)) +
        geom_point(aes(col = Water))+
        theme(text = element_text(size=12),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
              axis.title = element_text(size = 12))+
        ylab(bquote(italic(P.~secunda)~Biomass~(BRTE-None))) +
        xlab("Trait platicity")+
        geom_smooth(method = "lm", colour="black", size=0.5)


# Graph them together
ggarrange(f1, f2, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
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

p.mat <- cor.mtest(survival_trait[,7:ncol(survival_trait)])# matrix of p-values

corrplot(cor(survival_trait[,7:ncol(survival_trait)]), type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.01, insig = "blank")

ggplot(survival_trait%>%pivot_longer(cols = Length:LDMC, names_to = "trait", values_to = "value"), aes(x = value, y = POSE_survival_stem_count))+
  geom_point()+
  facet_wrap(~trait, scale = "free")
