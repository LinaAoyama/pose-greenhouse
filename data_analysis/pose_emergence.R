##Load R packages
library(tidyverse)
library(ggplot2)
library(SingleCaseES) 
library(ARPobservation)
install.packages("ARPobservation")

##Import data
source("data_compiling/compile_demography.R")
stems <- left_join(stemcount_dat, potID) #combine potID and stemcount data in one dataframe
leafarea<-left_join(leaftraits,potID)%>% left_join(.,biomass)%>% left_join(.,height)%>%
  mutate(SLA=FreshLeafArea_cm2/DryLeafWeight_g,LDMC=DryLeafWeight_g/FreshLeafWeight_g)

summarizeleaf<-leafarea%>%group_by(Population, Species)%>%filter(!is.na(Dry_Biomass_Weight_g), !is.na(SLA), !is.na(LDMC), !is.na(Height_cm))%>%summarise(meanheight=mean(Height_cm),seHeight=calcSE(Height_cm),meanSLA=mean(SLA),seSLA=calcSE(SLA),meanLDMC=mean(LDMC),seLDMC=calcSE(LDMC),meanDryBiomassWeight=mean(Dry_Biomass_Weight_g),seDryBiomassWeight=calcSE(Dry_Biomass_Weight_g))
summarizeseeds<-seeds%>%group_by(Population,Species)%>%filter(!is.na(Weight_g))%>%summarise(meanseed=mean(Weight_g),seseed=calcSE(Weight_g))
alltraits<-left_join(summarizeleaf,summarizeseeds)


##Visualize data
#scatter plots hypothesis #2 difference in traits by population
#Height
A_Height<- ggplot(summarizeleaf, aes(x = Population, y = meanheight)) + 
  geom_errorbar(aes(ymin=meanheight-seHeight,ymax=meanheight+seHeight))+
  #geom_jitter()+
  theme_classic()

#DryBiomassWeight
B_DryBiomassWeight<- ggplot(summarizeleaf, aes(x = Population, y = meanDryBiomassWeight)) + 
  geom_errorbar(aes(ymin=meanDryBiomassWeight-seDryBiomassWeight,ymax=meanDryBiomassWeight+seDryBiomassWeight))+
  #geom_jitter()+
  theme_classic()
#LDMC
C_LDMC<- ggplot(summarizeleaf, aes(x = Population, y = meanLDMC)) + 
  geom_errorbar(aes(ymin=meanLDMC-seLDMC,ymax=meanLDMC+seLDMC))+
  #geom_jitter()+
  theme_classic()
#SLA
D_SLA<- ggplot(summarizeleaf, aes(x = Population, y = meanSLA)) + 
  geom_errorbar(aes(ymin=meanSLA-seSLA,ymax=meanSLA+seSLA))+
  #geom_jitter()+
  theme_classic()
#Seeds
E_Seeds<- ggplot(summarizeseeds, aes(x = Population, y = meanseed)) + 
  geom_errorbar(aes(ymin=meanseed-seseed,ymax=meanseed+seseed))+
  #geom_jitter()+
  theme_classic()
   
install.packages("ggpubr")
  library(ggpubr)
ggarrange(A_Height, B_DryBiomassWeight, C_LDMC, D_SLA, E_Seeds, ncol=2, nrow=3, labels=c("Height","DryBiomassWeight","LDMC","SLA", "SeedWeight"))
  
#POSE_survival_stem_counts **changed col = water TO: col = Water**
ggplot(stems, aes(x = Water, y = POSE_survival_stem_counts)) + #x=Water, col=Water
  geom_point()+
  geom_jitter()+
  theme_bw()
geom_boxplot

#POSE_emergence_stem_counts
#boxplot
ggplot(stems, aes(x = Population, y = POSE_emergence_stem_counts)) +
  geom_boxplot()+
  geom_jitter()+
  theme_classic()
facet_wrap(~Water)

#POSE_survival_stem_counts boxplot
ggplot(stems, aes(x = Water, y = POSE_survival_stem_counts)) + #x=Water, col=Water
  geom_boxplot()+
  geom_jitter()+
  theme_classic()
facet_wrap(~Water)

#hypothesis three Does poa resistance to b tectorum change with water availability/competition?
#POSE_survival_stem_Counts boxplot with competition
ggplot(stems, aes(x = Competition, y = POSE_survival_stem_counts)) + #x=Water, col=Water
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  #facet_grid(~Population)
  facet_grid(vars(Water), vars(Population)) #this figure to show across population
  #facet_wrap(~Water)+ # this figure to show wet vs dry
  #facet_wrap(~Population)

#to show wet and dry -- show log response ratio
#stem counts and biomass --  %>% group_by(Population, Water, Competition)
#with(stems %>% filter(POSE_survival_stem_counts != "NA"),
#logRespRatio(observations = POSE_survival_stem_counts, phase = Water, base_level = "Dry", conf_level = 0.95,
#             bias_correct = TRUE, exponentiate = FALSE))

#lm(Y~X)  for NONE
lm_out <- lm(POSE_survival_stem_counts ~ Water,data=stems %>%filter(Competition == "None"))
beta_0 <- coef(lm_out)[1] 
beta_1 <- coef(lm_out)[2] 
#NONE object log response ratio:
lrr_None <- log(beta_0 + beta_1) - log(beta_0)

#lm(Y~X) for BRTE
lm_out <- lm(POSE_survival_stem_counts ~ Water,data=stems %>%filter(Competition == "BRTE"))
beta_0 <- coef(lm_out)[1] 
beta_1 <- coef(lm_out)[2] 
#BRTE object log response ratio:
lrr_BRTE <- log(beta_0 + beta_1) - log(beta_0)

#combine objects of BRTE/NONE log rr
cbind(lrr_None, lrr_BRTE)

logresponseratio <-cbind(lrr_None, lrr_BRTE)

#plot

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

  ##figure 2 with facets as dry wet, with competition
  ggplot(stems, aes(x = Population, y = POSE_emergence_stem_counts, fill=Competition)) + #col=water ##OR fill = Water##
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0))+ #rotate x axis labels
    #geom_point()+
    #geom_jitter()+
    #theme_classic()+
    facet_grid(~Water) 
  
  ggplot(stems, aes(x = Population, y = POSE_emergence_stem_counts, fill=Water)) + #col=water ##OR fill = Water##
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0))+ #rotate x axis labels
    #geom_point()+
    #geom_jitter()+
    #theme_classic()+
    facet_grid(~Competition) 
#
#y=biomass, x = populations ## this one with wet/dry/competition by using facet_Grid
#
  ggplot(leafarea, aes(x = Population, y = Dry_Biomass_Weight_g, fill=Competition)) + #col=water ##OR fill = Water##
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0))+ #rotate x axis labels
    #geom_point()+
    #geom_jitter()+
    #theme_classic()+
    facet_grid(~Water) 
  
  ggplot(leafarea, aes(x = Population, y = Dry_Biomass_Weight_g, fill=Water)) + #col=water ##OR fill = Water##
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0))+ #rotate x axis labels
    #geom_point()+
    #geom_jitter()+
    #theme_classic()+
    facet_grid(~Competition) 
  
  
  
#biomass by competition -- this is hard to interpret...
ggplot(leafarea, aes(x = Competition, y = Dry_Biomass_Weight_g)) + #x=Water, col=Water
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  #facet_grid(~Population)
  facet_grid(vars(Water), vars(Population)) #added population/competition
  
#seed weight #this needs work...
  ggplot(seeds, aes(x = Population, y = Weight_g))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0))
    #geom_point()+
    geom_jitter()+
    theme_classic()
  
  #height #no differences in hieght across population or water treatment?
  ggplot(leafarea, aes(x = Population, y = Height_cm, fill=Water))+
    geom_boxplot()+
    #geom_point()+
    #geom_jitter()+
    theme_classic()
  #facet_grid(~Water) 
  facet_grid(vars(Water), vars(Population)) #added population/competition
  
#height across population/water/competition  
ggplot(leafarea, aes(x = Competition, y = Height_cm)) + #x=Water, col=Water
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  #facet_grid(~Population)
  facet_grid(vars(Water), vars(Population)) #added population/competition
  
#sla, height, seed mass, ldmc


#figure 4 plasticity coefficient of variation--- traits from leafarea object:  SLA, LDMC
 traits_cv<-leafarea %>%
  filter(!is.na(Dry_Biomass_Weight_g), !is.na(SLA), !is.na(LDMC), !is.na(Height_cm))%>%
  group_by(Population) %>%
  summarise(dry_biomass_weight=sd(Dry_Biomass_Weight_g)/mean(Dry_Biomass_Weight_g)*100,
            SLA=sd(SLA)/mean(SLA)*100,
            LDMC=sd(LDMC)/mean(LDMC)*100,
            height=sd(Height_cm)/mean(Height_cm)*100)%>%
  gather(key="traits", value="cv",-Population)
  
  #plot ## This runs error "couldnt find function "scale_col_manual""?
ggplot(traits_cv, aes(x = traits, y = cv, col=Population)) + #col=water
    #geom_boxplot()+
  #geom_point(size=2)+
  geom_jitter()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0))+ #rotate x axis labels
 # theme_classic()
  #facet_wrap(~Water)
 scale_color_manual(name = "", values = c("#ff1f5b", "#00cd6c", "#009ade", "#af58ba", "#ffc61e", "#f28522"))
##b33000", "#ff4500",  "#ED7D31", "#5B9BD5", "#FF5733", "#70A62F"
###634071", "#303472",  "#395c45", "#7c7948", "#7d5039", "#6f2a2d
##brown:844d36", "dark grey:#474853",  "light blue;#86b3d1", "light grey:#aaa0a0", "tan:#8e8268", "green:#425a07
