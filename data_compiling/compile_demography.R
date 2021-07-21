library(tidyverse)
potID <- read_csv(paste(datpath, "/Greenhouse_potID.csv", sep = "")) #pot info and treatment
stemcount <- read_csv(paste(datpath, "/Greenhouse_seedling_stem_counts.csv", sep="")) #seedling counts
rametcount <- read_csv(datpath,"/Greenhouse_adult_ramet_counts.csv", sep="") #adult ramet counts
leaftraits <- read_csv(paste(datpath,"/Greenhouse_leaf_traits.csv", sep="")) #seedling fresh and dry leaf weight and area
biomass <- read_csv(paste(datpath,"/Greenhouse_aboveground_biomass.csv", sep="")) #aboveground BRTE and POSE dry weight
height <- read_csv(paste(datpath,"/Greenhouse_height.csv", sep="")) #seedling height
root <- read_csv(paste(datpath,"Greenhouse_roots.csv", sep = "")) #seedling root length, avg diameter, etc.
#seeds<- read_csv(paste(datpath,"/Greenhouse_POSE_seed_weight.csv", sep=""))