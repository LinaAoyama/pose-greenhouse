##Load R packages
library(tidyverse)
library(ggplot2)

##Import data
source("data_compiling/compile_demography.R")
stems <- left_join(stemcount_dat, potID) #combine potID and stemcount data in one dataframe
leafarea<-left_join(leaftraits,potID)%>% left_join(.,biomass)%>%
  mutate(SLA=FreshLeafArea_cm2/DryLeafWeight_g,LDMC=DryLeafWeight_g/FreshLeafWeight_g)


##Visualize data
#scatter plots 
#POSE_emergence_stem_counts
ggplot(stems, aes(x = Population, y = POSE_emergence_stem_counts)) +
  geom_point()+
  geom_jitter()+
  theme_classic()
  facet_wrap(~Water)

#POSE_survival_stem_counts **changed col = water TO: col = Water**
ggplot(stems, aes(x = Water, y = POSE_survival_stem_counts)) + #x=Water, col=Water
  geom_point()+
  geom_jitter()+
  theme_bw()
geom_boxplot

#POSE_emergence_stem_counts
ggplot(stems, aes(x = Population, y = POSE_emergence_stem_counts)) +
  geom_point()+
  geom_jitter()+
  theme_classic()
facet_wrap(~Water)
#can you think of other kinds of plots?

#formula to calculate standard error
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#calculate mean and se of stem counts per population
summary_stems <- stems %>% #object
  filter(POSE_survival_stem_counts != "NA") %>%
  group_by(Population, Water) %>%
  summarize(mean_POSE_em = mean(POSE_emergence_stem_counts), #mean(columnname)
            se_POSE_em = calcSE(POSE_emergence_stem_counts),
            mean_POSE_su = mean(POSE_survival_stem_counts),
            se_POSE_su = calcSE(POSE_survival_stem_counts),
            mean_BRTE = mean(BRTE_stem_counts),
            se_BRTE = calcSE(BRTE_stem_counts))



#bar plots
ggplot(summary_stems, aes(x = Population, y = mean_POSE_em))+ 
  #geom_
  geom_point()+
  geom_errorbar(aes(ymin = mean_POSE_em-se_POSE_em, ymax = mean_POSE_em+se_POSE_em), width = 0.4, alpha = 0.9, size = 1)+
  theme_bw()

#
#figure 2
#y= poa_growth (stems?), x = populations
#
ggplot(stems, aes(x = Population, y = POSE_emergence_stem_counts, fill=Water)) + #col=water
  geom_boxplot()+
  #geom_point()+
  #geom_jitter()+
  theme_classic()
  #facet_wrap(~Water)
  facet_grid(~Water) #how do I make the boxplot separate by wet/dry? this doesn't seem to work...
  
#OR#

#
#y=biomass, x = populations
#
 
#figure 4 plasticity coefficient of variation--- traits from leafarea object: seed weight, SLA, LDMC
 traits_cv<-leafarea %>%
  #group_by(Population) %>%
  #summarise(new column name cv_roots=sd(roots)/mean(roots))*100,cv_seedweight=sd(seedweight...)
  #
  
  #%>%gather(key="traits", value="cv",-Population)
  
  #plot
  ggplot(traits_cv, aes(x = traits, y = cv, col=Population)) + #col=water
    #geom_boxplot()+
    geom_point()+
    #geom_jitter()+
    theme_classic()
  #facet_wrap(~Water)
  
 #seeds plot separate

