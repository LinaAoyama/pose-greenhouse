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
summary(aov(BRTE~Life_stage, data = biomass_dat %>% filter(Competition == "BRTE")))
summary(aov(BRTE_stem_count~Life_stage, data = demography %>% filter(Competition == "BRTE")))

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
                    ylab(bquote(italic(B.tectorum)~Biomass~(g)))+
                    annotate("text", label = c("***"), y = 0.55, x = 1.5, size = 18)

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
                  ylab(bquote(italic(B.tectorum)~Germination~Rate))+
                  annotate("text", label = c("NA"), y = 0.95, x = 1.5, size = 10)

ggarrange( fig_biomass, fig_counts, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
           font.label = list(size = 18))

#-----------------------------------------------------------#
# 2. Are there any differences in survival rates of seedling POSE by population?

# Reorder the populations by wet to dry sites
growth$Population <- ordered(as.factor(growth$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                      "Water Canyon",  "Reno"))
biomass_dat$Population <- ordered(as.factor(biomass_dat$Population), levels = c("Butte Valley","Steens","EOARC", "Gund",
                                                                      "Water Canyon",  "Reno"))
# Calculate mean POSE stem counts by population
mu <- growth %>%
  filter(Water == "Wet") %>%
  group_by(Population, Competition) %>%
  summarise(mean = mean(POSE_survival_stem_count, na.rm = TRUE))

# Calculate survival rates
summary_seedling <- growth %>%
  drop_na()%>%
  filter(Water == "Wet") %>%
  group_by(Population, Competition) %>%
  summarise(mean = mean(POSE_survival_stem_count/25),
              se = se(POSE_survival_stem_count/25))

# Calculate mean BRTE biomass
summary_BRTE <- biomass_dat %>%
  drop_na()%>%
  filter(Competition == "BRTE")%>%
  filter(Life_stage == "seedling") %>%
  group_by(Population, Water) %>%
  summarise(mean = mean(BRTE),
            se = se(BRTE))
  
# Distribution of POSE seedling counts by population and BRTE competition - EOARC, Steens, and Water Canyon resisted BRTE
fig_density <- ggplot(growth%>%filter(Water == "Wet"), aes(x = POSE_survival_stem_count, fill = Competition)) +
                    geom_density(alpha = 0.2) +
                    theme_bw(base_size = 15)+
                    facet_grid(~Population)+
                    #geom_histogram(aes(y=..density..), alpha = 0.4, position = "identity")+
                    ylab(bquote(Density))+
                    xlab(bquote(italic(P.secunda)~stem~count))+
                    geom_vline(data = mu, aes(xintercept = mean, color = Competition))

# ggplot(growth%>%filter(Water == "Wet"), aes(x = POSE_survival_stem_count/25, fill = Competition)) +
#   geom_density(alpha = 0.2) +
#   theme_bw(base_size = 15)+
#   facet_grid(~Population)+
#   #geom_histogram(aes(y=..density..), alpha = 0.4, position = "identity")+
#   ylab(bquote(Density))+
#   xlab(bquote(italic(P.secunda)~Establishment~Rate))

# seeding POSE survival rate by population - EOARC, Gund, Water Canyon, and Reno more resistant to BRTE than Butte Valley and Steens
fig_establish <- ggplot(summary_seedling, aes(x = Population, y = mean, col = Competition)) +
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.position = c(0.8, 0.8), 
                        axis.title = element_text(size = 15))+
                  geom_point(position = position_dodge(width = 0.5))+
                  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                  ylab(bquote(italic(P.~secunda)~Survival~Rate)) +
                  annotate("text", label = c("A", "B", "B", "C", "B", "B"), y = 0.6, x = c(1, 2, 3, 4, 5, 6))

# Stats for POSE survival rate by pop
TukeyHSD(aov(POSE_survival_stem_count~Population, data = growth%>%
              filter(Water == "Wet")))

# Check stats for Butte Valley and Steens -- difference between BRTE and None not statistically significant
TukeyHSD(aov(POSE_survival_stem_count ~ Competition, data = growth%>%filter(Water == "Wet")%>%
               filter(Population == "Steens")))
TukeyHSD(aov(POSE_survival_stem_count ~ Competition, data = growth%>%filter(Water == "Wet")%>%
               filter(Population == "Butte Valley")))

# BRTE biomass by population 
fig_brte <- ggplot(summary_BRTE%>%filter(Water == "Wet"), aes(x = Population, y = mean))+
  geom_point() +
  theme_classic(base_size = 15) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2)+
  ylab(bquote(italic(B.tectorum)~Biomass~(g)))

# Stats for BRTE biomass by pop
TukeyHSD(aov(BRTE~Population, data = biomass_dat%>%
              filter(Life_stage == "seedling")%>%filter(Water == "Wet")%>%filter(Competition == "BRTE")))

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
                      legend.position = "top", 
                      axis.title = element_text(size = 15))+
                ylab(bquote(italic(P.secunda)~LRR~Survival~Rate))+
                geom_hline(yintercept = 0, linetype = "dashed")+
                scale_color_manual(name = "Water Treatment", values = c( "#E69F00", "#999999"))+
                annotate("text", label = c("*", "***", "*"), x = c(2, 5, 6), y = 1.1, size = 10)

# BRTE biomass by population and water availability 
fig_brte_water <- ggplot(summary_BRTE, aes(x = Population, y = mean, col = Water))+
                      geom_point(position = position_dodge(width = 0.5))+
                      geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1, position = position_dodge(width = 0.5))+
                      theme(text = element_text(size=15),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            legend.position = c(.7,0.9), 
                            axis.title = element_text(size = 15))+
                      ylab(bquote(italic(B.tectorum)~Biomass))+
                      scale_color_manual(name = "Water Treatment", values = c( "#E69F00", "#999999"))

# Graph them together
ggarrange(fig_LRR, fig_brte_water, ncol = 1, nrow = 2, labels = c("(a)", "(b)"),
          font.label = list(size = 15), legend = "top", heights = c(2, 2))
