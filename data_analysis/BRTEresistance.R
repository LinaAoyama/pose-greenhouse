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
growth <- inner_join(potID, stemcount) %>% filter(Population != "Gund")
ramets <- inner_join(potID, rametcount) %>% filter(Population != "Gund")
demography <- full_join(growth, ramets) %>% filter(Population != "Gund")
biomass_dat <- inner_join(potID, biomass) %>%
  filter(Population != "Gund") %>%
  pivot_wider(names_from = Species, values_from = Dry_Biomass_Weight_g) %>%
  replace(is.na(.),0)

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
# 2. Are there any differences in BRTE resistance of seedling POSE by population?

# Reorder the populations by wet to dry sites
growth$Population <- ordered(as.factor(growth$Population), levels = c("Butte Valley","Steens","EOARC", 
                                                                      "Water Canyon",  "Reno"))
growth$Water <- ordered(as.factor(growth$Water), levels = c("Wet", "Dry"))
growth$Treatment <- apply(growth[ ,3:4 ] , 1 , paste , collapse = "-" )
biomass_dat$Population <- ordered(as.factor(biomass_dat$Population), levels = c("Butte Valley","Steens","EOARC", 
                                                                      "Water Canyon",  "Reno"))
biomass_dat$Water <- ordered(as.factor(biomass_dat$Water), levels = c("Wet", "Dry"))
biomass_dat$Treatment <- apply(biomass_dat[ ,3:4 ] , 1 , paste , collapse = "-" )
# Calculate survival rates
summary_seedling <- growth %>%
  group_by(Population, Treatment) %>%
  summarise(mean = mean(POSE_survival_stem_count/25),
              se = se(POSE_survival_stem_count/25))

# Calculate mean biomass
summary_biomass <- biomass_dat %>%
  filter(Life_stage == "seedling") %>%
  dplyr::select(PotID, Population, Treatment, POSE) %>%
  drop_na()%>%
  group_by(Population, Treatment) %>%
  summarise(mean = mean(POSE), se = se(POSE))
  
# seedling POSE survival rate by population 
fig_establish <- ggplot(summary_seedling, aes(x = Treatment, y = mean, col = Population))+
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.position = c(0.2, 0.8), 
                        axis.title = element_text(size = 15))+
                  geom_point(position = position_dodge(width = 0.5))+
                  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                  ylab(bquote(italic(P.~secunda)~Establishment)) 

# Stats for POSE survival rate 
summary(aov(POSE_survival_stem_count ~ Population*Treatment, data = growth)) #both pop and treatment differences are significant
summary(lme(POSE_survival_stem_count ~ Treatment*Population, random = ~ 1|Replicate, data = growth))

# seedling POSE biomass by population
fig_biomass <- ggplot(summary_biomass, aes(x = Treatment, y = mean, col = Population)) +
                        theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.position = c(0.2, 0.8), 
                        axis.title = element_text(size = 15))+
                  geom_point(position = position_dodge(width = 0.5))+
                  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                  scale_y_log10()+
                  ylab(bquote(italic(P.~secunda)~Biomass~(g))) 

# Stats for interaction of pop and treatment on POSE biomass 
summary(aov(POSE ~ Treatment*Population, data = biomass_dat%>%filter(Life_stage == "seedling")))  #sig treatment differences but no pop or interaction
summary(lme(POSE ~ Treatment*Population, random = ~ 1|Replicate, data = biomass_dat%>%filter(Life_stage == "seedling")))

# Stats for pop differences
TukeyHSD(aov(POSE ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")%>%filter(Treatment == "BRTE-Dry")))
TukeyHSD(aov(POSE ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")%>%filter(Treatment == "BRTE-Wet")))
TukeyHSD(aov(POSE ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")%>%filter(Treatment == "None-Dry")))
TukeyHSD(aov(POSE ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")%>%filter(Treatment == "None-Wet")))

# Graph them together
ggarrange( fig_establish, fig_biomass, ncol = 1, nrow = 2, labels = c("(a)", "(b)"),
           font.label = list(size = 15), legend = "right", common.legend = TRUE)


# Distribution of POSE seedling counts by population and BRTE competition - EOARC, Steens, and Water Canyon resisted BRTE
# # Calculate mean POSE stem counts by population
# mu <- growth %>%
#   group_by(Population, Competition, Water) %>%
#   summarise(mean = mean(POSE_survival_stem_count, na.rm = TRUE))
# ggplot(growth%>%filter(Water == "Wet"), aes(x = POSE_survival_stem_count, fill = Competition)) +
#                     geom_density(alpha = 0.2) +
#                     theme_bw(base_size = 15)+
#                     facet_grid(~Population)+
#                     #geom_histogram(aes(y=..density..), alpha = 0.4, position = "identity")+
#                     ylab(bquote(Density))+
#                     xlab(bquote(italic(P.secunda)~stem~count))+
#                     geom_vline(data = mu, aes(xintercept = mean, color = Competition))

# ggplot(growth%>%filter(Water == "Wet"), aes(x = POSE_survival_stem_count/25, fill = Competition)) +
#   geom_density(alpha = 0.2) +
#   theme_bw(base_size = 15)+
#   facet_grid(~Population)+
#   #geom_histogram(aes(y=..density..), alpha = 0.4, position = "identity")+
#   ylab(bquote(Density))+
#   xlab(bquote(italic(P.secunda)~Establishment~Rate))


#BRTE biomass as metric of BRTE resistance?
# 
# # Calculate mean BRTE biomass
# summary_BRTE <- biomass_dat %>%
#   drop_na()%>%
#   filter(Competition == "BRTE")%>%
#   filter(Life_stage == "seedling") %>%
#   group_by(Population, Water) %>%
#   summarise(mean = mean(BRTE),
#             se = se(BRTE))
# 
# # BRTE biomass by population 
# fig_brte <- ggplot(summary_BRTE%>%filter(Water == "Wet"), aes(x = Population, y = mean))+
#   geom_point() +
#   theme_classic(base_size = 15) +
#   geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2)+
#   ylab(bquote(italic(B.tectorum)~Biomass~(g)))
# 
# # Stats for BRTE biomass by pop
# TukeyHSD(aov(BRTE~Population, data = biomass_dat%>%
#               filter(Life_stage == "seedling")%>%filter(Water == "Wet")%>%filter(Competition == "BRTE")))

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
