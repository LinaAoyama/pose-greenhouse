# Preliminary data exploration
# Cheatgrass (BRTE) resistance by Sandberg bluegrass (POSE) by population and by water availability

# Set pathway first

# Data
source("data_compiling/compile_demography.R")

# Packages
library(tidyverse)
library(ggplot2)
library(ggpubr)

# Prep data
growth <- inner_join(potID, stemcount)
ramets <- inner_join(potID, rametcount)
demography <- full_join(growth, ramets)
biomass_dat <- inner_join(potID, biomass) %>%
  pivot_wider(names_from = Species, values_from = Dry_Biomass_Weight_g)

#----------------------------------------------------------#
# 1. Are adult POSE better at resisting BRTE than seedling POSE?

# Adult vs seedling POSE - adult POSE are better than seedling POSE at limiting BRTE growth but not germination
fig_biomass <- ggplot(biomass_dat %>% filter(Competition == "BRTE"), aes(y = BRTE, x = Life_stage)) +
                    geom_boxplot()+
                    theme(text = element_text(size=20),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          legend.position = "none", 
                          axis.title = element_text(size = 18),
                          axis.title.x = element_blank())+
                    ylab(bquote(italic(B.tectorum)~Biomass~(g)))
fig_counts <- ggplot(demography %>% filter(Competition == "BRTE"), aes(y = BRTE_stem_count/50, x = Life_stage))+
                  geom_boxplot()+
                  theme(text = element_text(size=20),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.position = "none", 
                        axis.title = element_text(size = 18),
                        axis.title.x = element_blank())+
                  ylab(bquote(italic(B.tectorum)~Germination~Rate))
ggarrange( fig_biomass, fig_counts, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
           font.label = list(size = 18))

#-----------------------------------------------------------#
# 2. Are there any differences in growth rates of POSE by population?

# Reorder the populations by wet to dry sites
ramets$Population <- ordered(as.factor(ramets$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                               "Water Canyon",  "Reno"))
growth$Population <- ordered(as.factor(growth$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                      "Water Canyon",  "Reno"))
# adult POSE growth rate by population
ggplot(ramets, aes(x = Population, y = POSE_ramet_count_post/POSE_ramet_count_pre, fill = Competition)) +
  geom_boxplot()

# seeding POSE growth rate by population
ggplot(growth%>%filter(Water == "Wet"), aes(x = Population, y = POSE_survival_stem_count/25, fill = Competition)) +
  geom_boxplot()

# Distribution of POSE seedling counts by population and BRTE competition - EOARC, Steens, and Water Canyon resisted BRTE
ggplot(growth, aes(x = POSE_survival_stem_count, fill = Competition)) +
  geom_density(alpha = 0.4) +
  theme_bw()+
  facet_wrap(~Population)+
  geom_histogram(aes(y=..density..), alpha = 0.4, position = "identity")



#-----------------------------------------------------------#
# 3. Does POSE's resistance to BRTE change by water availability?


