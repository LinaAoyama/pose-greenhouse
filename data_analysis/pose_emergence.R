##Load R packages
library(tidyverse)
library(ggplot2)

##Import data
source("data_compiling/compile_demography.R")
stems <- left_join(stemcount_dat, potID) #combine potID and stemcount data in one dataframe

##Visualize data
#scatter plots
ggplot(stems, aes(x = Population, y = POSE_emergence_stem_counts)) +
  geom_point()+
  geom_jitter()+
  theme_bw()
#  facet_wrap(~Water)

ggplot(stems, aes(x = Population, y = POSE_survival_stem_counts)) +
  geom_point()+
  geom_jitter()+
  theme_bw()

#can you think of other kinds of plots?

#formula to calculate standard error
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#calculate mean and se of stem counts per population
summary_stems <- stems %>%
  filter(POSE_survival_stem_counts != "NA") %>%
  group_by(Population) %>%
  summarize(mean_POSE_em = mean(POSE_emergence_stem_counts),
            se_POSE_em = calcSE(POSE_emergence_stem_counts),
            mean_POSE_su = mean(POSE_survival_stem_counts),
            se_POSE_su = calcSE(POSE_survival_stem_counts),
            mean_BRTE = mean(BRTE_stem_counts),
            se_BRTE = calcSE(BRTE_stem_counts))

#bar plots
ggplot(summary_stems, aes(x = Population, y = mean_POSE_em))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_POSE_em-se_POSE_em, ymax = mean_POSE_em+se_POSE_em), width = 0.4, alpha = 0.9, size = 1)+
  theme_bw()
