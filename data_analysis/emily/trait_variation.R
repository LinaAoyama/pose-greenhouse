##Load packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
##Import data
source("data_compiling/compile_demography.R")
adultstem <- read_csv(paste(datpath, "/Greenhouse_adult_ramet_counts.csv", sep = ""))

##Combine data
leaf<-leaftraits %>%
  mutate(SLA=FreshLeafArea_cm2/DryLeafWeight_g,LDMC=DryLeafWeight_g/FreshLeafWeight_g) %>%
  filter(!is.na(SLA))
ht <- height %>%
  group_by(PotID) %>%
  summarise(height_cm = mean(Height_cm))
data <- left_join(potID, stemcount_dat) %>% left_join(., biomass %>% filter(Species == "POSE")) %>% left_join(., leaf) %>% left_join(., ht)
data$Population <- ordered(data$Population, levels = c("Butte Valley",  "Steens", "EOARC",  "Gund", "Water Canyon", "Reno"))

##Calculate growth rate 
#GR = ln((output+1)/(input+1))
adult <- adultstem %>%
  mutate(GR = log((POSE_ramet_count_post+1)/(POSE_ramet_count_pre+1))) %>%
  left_join(., potID, by = "PotID")

##Log response ratios of adult growth rate
lnRRadult <- adult %>%
  group_by(Competition, Water, Population) %>%
  summarise(meanGR = mean(GR, na.rm = T)) %>%
  pivot_wider(names_from = "Competition", values_from = "meanGR" ) %>%
  mutate(lnRR = log(BRTE/None))
lnRRadult$Population <- ordered(lnRRadult$Population, levels = c("Butte Valley",  "Steens", "EOARC",  "Gund", "Water Canyon", "Reno"))

##Calculate survival rate
seedling <- stemcount_dat %>%
  mutate(survivalrate = POSE_survival_stem_counts/25) %>%
  left_join(., potID, by = "PotID") 
  
##Log response ratios of survival rate
lnRRseedling <- seedling %>%
  group_by(Competition, Water, Population) %>%
  summarise(meanSurvival = mean(survivalrate, na.rm = T)) %>%
  pivot_wider(names_from = "Competition", values_from = "meanSurvival") %>%
  mutate(lnRR = log(BRTE/None))
lnRRseedling$Population <- ordered(lnRRseedling$Population, levels = c("Butte Valley",  "Steens", "EOARC",  "Gund", "Water Canyon", "Reno"))

##Visualize lnRR adult growth rate and seedling survival rate
a <- ggplot(lnRRadult, aes(x = Population, y = lnRR))+
  geom_bar(stat = "identity", fill = "#00BFC4") +
  theme_classic()+
  geom_hline(yintercept = 0) +
  labs(y = "lnRR ramet growth rate")
b <- ggplot(lnRRseedling, aes(x = Population, y = lnRR, fill = Water)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  geom_hline(yintercept = 0) +
  labs(y = "lnRR survival rate")

ggarrange(a, b)
##Visualize POSE seedling stem count by watering trt

ggplot(data %>% filter(Life_stage == "seedling"), aes(x = Population, y = POSE_survival_stem_counts, fill = Water)) +
  geom_boxplot()

##Visualize POSE seedling stem count by competition trt

ggplot(data %>% filter(Life_stage == "seedling"), aes(x = Population, y = POSE_survival_stem_counts, fill = Competition)) +
  geom_boxplot()

##Visualize POSE aboveground biomass by competition and watering trt

ggplot(data %>% filter(Life_stage == "seedling"), aes(x = Population, y = Dry_Biomass_Weight_g, fill = Competition)) +
  geom_boxplot()+
  facet_grid(~Water)

##Visualize POSE height by competition and watering trt

ggplot(data %>% filter(Life_stage == "seedling"), aes(x = Population, y = height_cm, fill = Competition)) +
  geom_boxplot()+
  facet_grid(~Water)

##Visualize POSE SLA by competition and watering trt

ggplot(data %>% filter(Life_stage == "seedling"), aes(x = Population, y = SLA, fill = Competition)) +
  geom_boxplot()+
  facet_grid(~Water)

##Visualize POSE LDMC by competition and watering trt

ggplot(data %>% filter(Life_stage == "seedling"), aes(x = Population, y = LDMC, fill = Competition)) +
  geom_boxplot()+
  facet_grid(~Water)




