# Preliminary data exploration
# Cheatgrass (BRTE) resistance by Sandberg bluegrass (POSE) by population and by water availability

# Set pathway first

# Data
source("data_compiling/compile_demography.R")

# Packages
library(tidyverse) #data wrangling
library(ggplot2) #plot
library(ggpubr) #combine plots
library(SingleCaseES) #log response ratios
library(nlme) #linear mixed effects
library(multcomp) #tukey

# Prep data
growth <- inner_join(potID, stemcount)
ramets <- inner_join(potID, rametcount)
demography <- full_join(growth, ramets)
biomass_dat <- inner_join(potID, biomass) %>%
  pivot_wider(names_from = Species, values_from = Dry_Biomass_Weight_g)

# Function for standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
} 

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
# 2. Are there any differences in establishment rates of seedling POSE by population?

# Reorder the populations by wet to dry sites
growth$Population <- ordered(as.factor(growth$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                      "Water Canyon",  "Reno"))

# Calculate mean
mu <- growth %>%
  filter(Water == "Wet") %>%
  group_by(Population, Competition) %>%
  summarise(mean = mean(POSE_survival_stem_count, na.rm = TRUE))

# Calculate establishment rates
summary_seedling <- growth %>%
  drop_na()%>%
  filter(Water == "Wet") %>%
  group_by(Population, Competition) %>%
  summarise(mean = mean(POSE_survival_stem_count/25),
              se = se(POSE_survival_stem_count/25))
  
# Distribution of POSE seedling counts by population and BRTE competition - EOARC, Steens, and Water Canyon resisted BRTE
fig_density <- ggplot(growth%>%filter(Water == "Wet"), aes(x = POSE_survival_stem_count, fill = Competition)) +
                    geom_density(alpha = 0.2) +
                    theme_bw(base_size = 15)+
                    facet_grid(~Population)+
                    #geom_histogram(aes(y=..density..), alpha = 0.4, position = "identity")+
                    ylab(bquote(Density))+
                    xlab(bquote(italic(P.secunda)~stem~count))+
                    geom_vline(data = mu, aes(xintercept = mean, color = Competition))

# seeding POSE establishment rate by population - EOARC, Gund, Water Canyon, and Reno more resistant to BRTE than Butte Valley and Steens
fig_establish <- ggplot(summary_seedling, aes(x = Population, y = mean, col = Competition)) +
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.position = c(0.8, 0.9), 
                        axis.title = element_text(size = 15))+
                  geom_point(position = position_dodge(width = 0.5))+
                  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                  ylab(bquote(italic(P.secunda)~Establishmeent~Rate))

# Check stats for Butte Valley and Steens -- difference between BRTE and None not statistically significant
TukeyHSD(aov(POSE_survival_stem_count ~ Competition, data = growth%>%filter(Water == "Wet")%>%
               filter(Population == "Steens")))
TukeyHSD(aov(POSE_survival_stem_count ~ Competition, data = growth%>%filter(Water == "Wet")%>%
               filter(Population == "Butte Valley")))

# Graph them together
ggarrange(fig_density, fig_establish, ncol = 1, nrow = 2, labels = c("(a)", "(b)"),
           font.label = list(size = 15), legend = "top", heights = c(1.5, 2))

#-----------------------------------------------------------#
# 3. Does POSE's resistance to BRTE change by water availability?

# Log-response ratios of BRTE competition treatment/control in seedling POSE
LRR_seedlings_wet <- rbind(LRRd(A_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Butte Valley")%>%filter(Competition == "None")%>%select(growth), 
                         B_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Butte Valley")%>%filter(Competition == "BRTE")%>%select(growth)),
                    LRRd(A_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Steens")%>%filter(Competition == "None")%>%select(growth), 
                         B_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Steens")%>%filter(Competition == "BRTE")%>%select(growth)),
                    LRRd(A_data = growth%>%filter(Water == "Wet")%>%filter(Population == "EOARC")%>%filter(Competition == "None")%>%select(growth), 
                         B_data = growth%>%filter(Water == "Wet")%>%filter(Population == "EOARC")%>%filter(Competition == "BRTE")%>%select(growth)),
                    LRRd(A_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Gund")%>%filter(Competition == "None")%>%select(growth), 
                         B_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Gund")%>%filter(Competition == "BRTE")%>%select(growth)),
                    LRRd(A_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Water Canyon")%>%filter(Competition == "None")%>%select(growth), 
                         B_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Water Canyon")%>%filter(Competition == "BRTE")%>%select(growth)),
                    LRRd(A_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Reno")%>%filter(Competition == "None")%>%select(growth), 
                         B_data = growth%>%filter(Water == "Wet")%>%filter(Population == "Reno")%>%filter(Competition == "BRTE")%>%select(growth))) %>%
              mutate(Population = c("Butte Valley", "Steens", "EOARC", 
                                    "Gund", "Water Canyon", "Reno"),
                     Water = "Wet")
LRR_seedlings_dry <- rbind(LRRd(A_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Butte Valley")%>%filter(Competition == "None")%>%select(growth), 
                                B_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Butte Valley")%>%filter(Competition == "BRTE")%>%select(growth)),
                           LRRd(A_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Steens")%>%filter(Competition == "None")%>%select(growth), 
                                B_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Steens")%>%filter(Competition == "BRTE")%>%select(growth)),
                           LRRd(A_data = growth%>%filter(Water == "Dry")%>%filter(Population == "EOARC")%>%filter(Competition == "None")%>%select(growth), 
                                B_data = growth%>%filter(Water == "Dry")%>%filter(Population == "EOARC")%>%filter(Competition == "BRTE")%>%select(growth)),
                           LRRd(A_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Gund")%>%filter(Competition == "None")%>%select(growth), 
                                B_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Gund")%>%filter(Competition == "BRTE")%>%select(growth)),
                           LRRd(A_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Water Canyon")%>%filter(Competition == "None")%>%select(growth), 
                                B_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Water Canyon")%>%filter(Competition == "BRTE")%>%select(growth)),
                           LRRd(A_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Reno")%>%filter(Competition == "None")%>%select(growth), 
                                B_data = growth%>%filter(Water == "Dry")%>%filter(Population == "Reno")%>%filter(Competition == "BRTE")%>%select(growth))) %>%
              mutate(Population = c("Butte Valley", "Steens", "EOARC", 
                                    "Gund", "Water Canyon", "Reno"),
                     Water = "Dry")
LRR_seedlings <- rbind(LRR_seedlings_wet, LRR_seedlings_dry)
LRR_seedlings$Population <- ordered(as.factor(LRR_seedlings$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                      "Water Canyon",  "Reno"))
fig_LRR <- ggplot(LRR_seedlings, aes(x = Population, y = Est, col = Water))+
                geom_point(position = position_dodge(width = 0.5))+
                geom_errorbar(aes(ymin = Est-SE, ymax = Est+SE), width = 0.2, alpha = 0.9, size = 1, position = position_dodge(width = 0.5))+
                theme(text = element_text(size=15),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position = c(.2,0.9), 
                      axis.title = element_text(size = 15))+
                ylab(bquote(log(BRTE/None)))+
                geom_hline(yintercept = 0, linetype = "dashed")+
                scale_color_manual(name = "Water Treatment", values = c( "#E69F00", "#999999"))

