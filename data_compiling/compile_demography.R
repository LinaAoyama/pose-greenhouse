library(tidyverse)
stemcount_dat <- read_csv(paste(datpath, "/Greenhouse_seedling_stem_counts.csv", sep="")) 
potID <- read_csv(paste(datpath, "/Greenhouse_potID.csv", sep = ""))
leaftraits <- read_csv(paste(datpath,"/Greenhouse_leaf_traits.csv", sep=""))
