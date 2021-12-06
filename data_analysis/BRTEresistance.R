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
  filter(Species == "POSE") %>%
  pivot_wider(names_from = Species, values_from = Dry_Biomass_Weight_g) %>%
  replace(is.na(.),0) %>%
  inner_join(., rootbiomass) %>%
  mutate(TotalBiomass = POSE + Root_Weight_g)

# Function for standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
} 

#----------------------------------------------------------#
# 1. Are adult POSE better at resisting BRTE than seedling POSE?

# Adult vs seedling POSE - adult POSE are better than seedling POSE at limiting BRTE growth but not germination
# summary(aov(BRTE~Life_stage, data = biomass_dat %>% filter(Competition == "BRTE")))
# summary(aov(BRTE_stem_count~Life_stage, data = demography %>% filter(Competition == "BRTE")))
# 
# fig_biomass <- ggplot(biomass_dat %>% filter(Competition == "BRTE"), aes(y = BRTE, x = Life_stage)) +
#                     geom_boxplot()+
#                     theme(text = element_text(size=20),
#                           panel.grid.major = element_blank(),
#                           panel.grid.minor = element_blank(),
#                           panel.background = element_blank(),
#                           axis.line = element_line(colour = "black"),
#                           legend.position = "none", 
#                           axis.title = element_text(size = 18),
#                           axis.title.x = element_blank())+
#                     ylab(bquote(italic(B.tectorum)~Biomass~(g)))+
#                     annotate("text", label = c("***"), y = 0.55, x = 1.5, size = 18)
# 
# fig_counts <- ggplot(demography %>% filter(Competition == "BRTE"), aes(y = BRTE_stem_count/50, x = Life_stage))+
#                   geom_boxplot()+
#                   theme(text = element_text(size=20),
#                         panel.grid.major = element_blank(),
#                         panel.grid.minor = element_blank(),
#                         panel.background = element_blank(),
#                         axis.line = element_line(colour = "black"),
#                         legend.position = "none", 
#                         axis.title = element_text(size = 18),
#                         axis.title.x = element_blank())+
#                   ylab(bquote(italic(B.tectorum)~Germination~Rate))+
#                   annotate("text", label = c("NA"), y = 0.95, x = 1.5, size = 10)
# 
# ggarrange( fig_biomass, fig_counts, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
#            font.label = list(size = 18))

#-----------------------------------------------------------#
# 2. Are there any differences in drought responses of seedling POSE by population?

# Reorder the populations by wet to dry sites
growth$Population <- ordered(as.factor(growth$Population), levels = c("Butte Valley","Steens","EOARC", 
                                                                      "Water Canyon",  "Reno"))
growth$Treatment <- apply(growth[ ,3:4 ] , 1 , paste , collapse = "-" )
growth$Treatment <- ordered(as.factor(growth$Treatment), levels = c("None-Wet", "None-Dry", "BRTE-Wet", "BRTE-Dry"))
biomass_dat$Population <- ordered(as.factor(biomass_dat$Population), levels = c("Butte Valley","Steens","EOARC", 
                                                                      "Water Canyon",  "Reno"))
biomass_dat$Treatment <- apply(biomass_dat[ ,3:4 ] , 1 , paste , collapse = "-" )
biomass_dat$Treatment <- ordered(as.factor(biomass_dat$Treatment), levels = c("None-Wet", "None-Dry", "BRTE-Wet", "BRTE-Dry"))

# Calculate mean survival rates by population
population_seedling <- growth %>%
  group_by(Population) %>%
  summarise(mean = mean(POSE_survival_stem_count/25),
            se = se(POSE_survival_stem_count/25))

fig_establish <-ggplot(population_seedling, aes(x = Population, y = mean))+
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.position = "none", 
                        axis.title = element_text(size = 15),
                        axis.title.x = element_blank())+
                  geom_point()+
                  ylab(bquote(Establishment~Rate))+
                  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                  annotate("text", label = c("*"), x = c(1), y = 0.45, size = 8)+
                  ylim(0, 0.5)


summary(aov(POSE_survival_stem_count ~ Population, data = growth%>%filter(Life_stage == "seedling")))
TukeyHSD(aov(POSE_survival_stem_count ~ Population, data = growth%>%filter(Life_stage == "seedling")))

# Calculate survival rates by treatment and population
trt_seedling <- growth %>%
  group_by(Treatment, Population) %>%
  summarise(mean = mean(POSE_survival_stem_count/25),
            se = se(POSE_survival_stem_count/25))

fig_establish_trt <-ggplot(trt_seedling, aes(x = Population, y = mean, col = Treatment))+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top", 
        axis.title = element_text(size = 15),
        axis.title.x = element_blank())+
  geom_point(position = position_dodge(width = 0.5))+
  ylab(bquote(Establishment~Rate))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, col = Treatment),width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
  annotate("text", label = c("*"), x = c(1), y = 0.45, size = 8)+
  scale_color_manual(values=c("#0240FF",  "#86BBE8",  "#FF8C07", "#FBD947"))

# # Calculate mean biomass by population
# population_biomass <- biomass_dat %>%
#   group_by(Population) %>%
#   summarise(mean = mean(TotalBiomass),
#             se = se(TotalBiomass))
# 
# fig_biomass <-ggplot(population_biomass, aes(x = Population, y = mean))+
#   theme(text = element_text(size=15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "none", 
#         axis.title = element_text(size = 15),
#         axis.title.x = element_blank())+
#   geom_point()+
#   ylab(bquote(Total~Biomass~(g)))+
#   geom_errorbar(aes(ymin = mean-se, ymax = mean+se),width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
#   ylim(0, 0.3)
# 
# summary(aov(TotalBiomass ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")))
# TukeyHSD(aov(TotalBiomass ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")))

# Calculate POSE's drought response by BRTE treatment
# Log-response ratios of survival rate
LRR_seedlings_None <- rbind(LRRd(A_data = growth%>%filter(Competition == "None")%>%filter(Population == "Butte Valley")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "None")%>%filter(Population == "Butte Valley")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count)),
                           LRRd(A_data = growth%>%filter(Competition == "None")%>%filter(Population == "Steens")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "None")%>%filter(Population == "Steens")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count)),
                           LRRd(A_data = growth%>%filter(Competition == "None")%>%filter(Population == "EOARC")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "None")%>%filter(Population == "EOARC")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count)),
                           LRRd(A_data = growth%>%filter(Competition == "None")%>%filter(Population == "Water Canyon")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "None")%>%filter(Population == "Water Canyon")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count)),
                           LRRd(A_data = growth%>%filter(Competition == "None")%>%filter(Population == "Reno")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "None")%>%filter(Population == "Reno")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count))) %>%
  mutate(Population = c("Butte Valley", "Steens", "EOARC", "Water Canyon", "Reno"),
         Competition = "None")
LRR_seedlings_BRTE <- rbind(LRRd(A_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "Butte Valley")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "Butte Valley")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count)),
                           LRRd(A_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "Steens")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "Steens")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count)),
                           LRRd(A_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "EOARC")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "EOARC")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count)),
                           LRRd(A_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "Water Canyon")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "Water Canyon")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count)),
                           LRRd(A_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "Reno")%>%filter(Water == "Wet")%>%dplyr::select(POSE_survival_stem_count), 
                                B_data = growth%>%filter(Competition == "BRTE")%>%filter(Population == "Reno")%>%filter(Water == "Dry")%>%dplyr::select(POSE_survival_stem_count))) %>%
  mutate(Population = c("Butte Valley", "Steens", "EOARC", 
                        "Water Canyon", "Reno"),
         Competition = "BRTE")
LRR_seedlings <- rbind(LRR_seedlings_None, LRR_seedlings_BRTE)
LRR_seedlings$Population <- ordered(as.factor(LRR_seedlings$Population), levels = c("Butte Valley","Steens","EOARC", 
                                                                                    "Water Canyon",  "Reno"))
fig_establish_DR <- ggplot(LRR_seedlings, aes(x = Population, y = Est, col = Competition))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = Est-SE, ymax = Est+SE), width = 0.2, alpha = 0.9, size = 1, position = position_dodge(width = 0.5))+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.5, 0.9), 
        legend.title = element_blank(),
        axis.title = element_text(size = 15))+
  ylab(bquote(LRR~Establishment~Rate))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(name = "Competition Treatment", values = c( "#777777", "#111111"))+
  annotate("text", label = c("*", "*"), x = c(2, 4), y = 0.5, size = 8)+
  ylim(-1.5, 0.9)

# # Log-response ratios of biomass
# LRR_biomass_None <- rbind(LRRd(A_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "Butte Valley")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "Butte Valley")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass)),
#                             LRRd(A_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "Steens")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "Steens")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass)),
#                             LRRd(A_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "EOARC")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "EOARC")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass)),
#                             LRRd(A_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "Water Canyon")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "Water Canyon")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass)),
#                             LRRd(A_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "Reno")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "None")%>%filter(Population == "Reno")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass))) %>%
#   mutate(Population = c("Butte Valley", "Steens", "EOARC", "Water Canyon", "Reno"),
#          Competition = "None")
# LRR_biomass_BRTE <- rbind(LRRd(A_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "Butte Valley")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "Butte Valley")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass)),
#                             LRRd(A_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "Steens")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "Steens")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass)),
#                             LRRd(A_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "EOARC")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "EOARC")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass)),
#                             LRRd(A_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "Water Canyon")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "Water Canyon")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass)),
#                             LRRd(A_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "Reno")%>%filter(Water == "Wet")%>%dplyr::select(TotalBiomass), 
#                                  B_data = biomass_dat%>%filter(Competition == "BRTE")%>%filter(Population == "Reno")%>%filter(Water == "Dry")%>%dplyr::select(TotalBiomass))) %>%
#   mutate(Population = c("Butte Valley", "Steens", "EOARC", 
#                         "Water Canyon", "Reno"),
#          Competition = "BRTE")
# LRR_biomass <- rbind(LRR_biomass_None, LRR_biomass_BRTE)
# LRR_biomass$Population <- ordered(as.factor(LRR_biomass$Population), levels = c("Butte Valley","Steens","EOARC", 
#                                                                                     "Water Canyon",  "Reno"))
# fig_biomass_DR <- ggplot(LRR_biomass, aes(x = Population, y = Est, col = Competition))+
#   geom_point(position = position_dodge(width = 0.5))+
#   geom_errorbar(aes(ymin = Est-SE, ymax = Est+SE), width = 0.2, alpha = 0.9, size = 1, position = position_dodge(width = 0.5))+
#   theme(text = element_text(size=15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "top", 
#         axis.title = element_text(size = 15))+
#   ylab(bquote(LRR~Total~Biomass))+
#   geom_hline(yintercept = 0, linetype = "dashed")+
#   scale_color_manual(name = "Competition Treatment", values = c( "#FF8C07", "#0240FF"))+
#   annotate("text", label = c("*", "*", "*", "*", "*"), x = c(1, 2, 3, 4, 5), y = 0.8, size = 8)+
#   ylim(-1.5, 0.9)

# Calculate survival rates by treatment and pop
summary_seedling <- growth %>%
  group_by(Population, Water, Competition, Treatment) %>%
  summarise(mean = mean(POSE_survival_stem_count/25),
              se = se(POSE_survival_stem_count/25))

# # Calculate mean total biomass by treatment and pop
# summary_biomass <- biomass_dat %>%
#   mutate(root_shoot_ratio = Root_Weight_g/POSE) %>%
#   dplyr::select(PotID, Population, Water, Competition, Treatment, TotalBiomass, root_shoot_ratio) %>%
#   drop_na()%>%
#   group_by(Population, Water, Competition, Treatment) %>%
#   summarise(mean = mean(TotalBiomass), se = se(TotalBiomass),
#             mean_ratio = mean(root_shoot_ratio), se_ratio = se(root_shoot_ratio))
# 
# seedling POSE survival rate by population
fig_establish_None <- ggplot(summary_seedling %>% filter(Competition == "None"), aes(x = Population, y = mean, col = Treatment))+
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.position = c(0.5, 0.8),
                        legend.title = element_blank(),
                        axis.title = element_text(size = 15),
                        axis.title.x = element_blank())+
                  geom_point(position = position_dodge(width = 0.5))+
                  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                  ylab(bquote(Establishment~Rate)) +
                  scale_color_manual(values=c( "#0240FF",  "#86BBE8"))+
                  annotate("text", label = c( "*", "*", "*"), x = c(2, 4, 5), y = 0.6, size = 8)
fig_establish_BRTE <- ggplot(summary_seedling %>% filter(Competition == "BRTE"), aes(x = Population, y = mean, col = Treatment))+
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.position = c(0.5, 0.8),
                        legend.title = element_blank(),
                        axis.title = element_text(size = 15),
                        axis.title.x = element_blank())+
                  geom_point(position = position_dodge(width = 0.5))+
                  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
                  ylab(bquote(Establishment~Rate)) +
                  scale_color_manual(values=c( "#FF8C07", "#FBD947"))+
                  annotate("text", label = c( "*"), x = c( 5), y = 0.6, size = 8)
# Stats for POSE survival rate
summary(aov(POSE_survival_stem_count ~ Population*Treatment, data = growth)) #both pop and treatment differences are significant
TukeyHSD(aov(POSE_survival_stem_count ~ Treatment, data = growth%>%filter(Life_stage == "seedling") %>%filter(Population == "Butte Valley")))
TukeyHSD(aov(POSE_survival_stem_count ~ Treatment, data = growth%>%filter(Life_stage == "seedling") %>%filter(Population == "Steens")))
TukeyHSD(aov(POSE_survival_stem_count ~ Treatment, data = growth%>%filter(Life_stage == "seedling") %>%filter(Population == "Water Canyon")))
TukeyHSD(aov(POSE_survival_stem_count ~ Treatment, data = growth%>%filter(Life_stage == "seedling") %>%filter(Population == "EOARC")))
TukeyHSD(aov(POSE_survival_stem_count ~ Treatment, data = growth%>%filter(Life_stage == "seedling") %>%filter(Population == "Reno")))
summary(lme(POSE_survival_stem_count ~ Treatment*Population, random = ~ 1|Replicate, data = growth))

# # seedling POSE total biomass by population
# fig_biomass_None <- ggplot(summary_biomass%>% filter(Competition == "None"), aes(x = Population, y = mean, col = Treatment)) +
#                         theme(text = element_text(size=15),
#                         panel.grid.major = element_blank(),
#                         panel.grid.minor = element_blank(),
#                         panel.background = element_blank(),
#                         axis.line = element_line(colour = "black"),
#                         legend.position = "none", 
#                         axis.title = element_text(size = 15))+
#                   geom_point(position = position_dodge(width = 0.5))+
#                   geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
#                   scale_y_log10()+
#                   ylab(bquote(Total~Biomass~(g))) +
#                   scale_color_manual(values=c("#56B4E9","#E69F00" ))
# fig_biomass_BRTE <- ggplot(summary_biomass%>% filter(Competition == "BRTE"), aes(x = Population, y = mean, col = Treatment)) +
#   theme(text = element_text(size=15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "none", 
#         axis.title = element_text(size = 15))+
#   geom_point(position = position_dodge(width = 0.5))+
#   geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
#   scale_y_log10()+
#   ylab(bquote(Total~Biomass~(g))) +
#   scale_color_manual(values=c("#6A0DAD", "#999999" ))

# Graph them together
ggarrange(fig_establish, fig_establish_None, fig_establish_BRTE,   fig_establish_DR,  
          ncol = 1, nrow = 4, labels = c("(a)", "(b)", "(c)", "(d)"),
          font.label = list(size = 15), align = "v", heights = c(1,1, 1, 1.2))

# 
# # Stats for interaction of pop and treatment on POSE biomass 
# summary(aov(POSE ~ Treatment*Population, data = biomass_dat%>%filter(Life_stage == "seedling")))  #sig treatment differences but no pop or interaction
# summary(lme(POSE ~ Treatment*Population, random = ~ 1|Replicate, data = biomass_dat%>%filter(Life_stage == "seedling")))
# 
# # Stats for pop differences
# TukeyHSD(aov(POSE ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")%>%filter(Treatment == "BRTE-Dry")))
# TukeyHSD(aov(POSE ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")%>%filter(Treatment == "BRTE-Wet")))
# TukeyHSD(aov(POSE ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")%>%filter(Treatment == "None-Dry")))
# TukeyHSD(aov(POSE ~ Population, data = biomass_dat%>%filter(Life_stage == "seedling")%>%filter(Treatment == "None-Wet")))
# 
# # seedling POSE root to shoot ratio by population
# fig_ratio <- ggplot(summary_biomass, aes(x = Population, y = mean_ratio, col = Treatment)) +
#   theme(text = element_text(size=15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = c(0.2, 0.8), 
#         axis.title = element_text(size = 15))+
#   geom_point(position = position_dodge(width = 0.5))+
#   geom_errorbar(aes(ymin = mean_ratio-se_ratio, ymax = mean_ratio+se_ratio), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
#   scale_y_log10()+
#   ylab(bquote(Root/Shoot~Ratio)) +
#   scale_color_manual(values=c("#56B4E9","#E69F00", "#6A0DAD", "#999999" ))

# # Seedling mortality
# summary_mortality <- growth %>%
#   mutate(delta =  POSE_survival_stem_count- POSE_emergence_stem_count) %>%
#   group_by(Population, Treatment) %>%
#   summarise(mean = mean(delta),
#             se = se(delta))
# 
# # seedling POSE mortality by population
# ggplot(summary_mortality, aes(x = Treatment, y = mean, col = Population)) +
#   theme(text = element_text(size=15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = c(0.2, 0.2), 
#         axis.title = element_text(size = 15))+
#   geom_point(position = position_dodge(width = 0.5))+
#   geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
#   geom_hline(yintercept=0)+
#   ylab(bquote(italic(P.~secunda)~Delta~("9wks-6wks")) )

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

# # BRTE biomass by population and water availability 
# fig_brte_water <- ggplot(summary_BRTE, aes(x = Population, y = mean, col = Water))+
#   geom_point(position = position_dodge(width = 0.5))+
#   geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1, position = position_dodge(width = 0.5))+
#   theme(text = element_text(size=15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = c(.7,0.9), 
#         axis.title = element_text(size = 15))+
#   ylab(bquote(italic(B.tectorum)~Biomass))+
#   scale_color_manual(name = "Water Treatment", values = c( "#E69F00", "#999999"))
# 
