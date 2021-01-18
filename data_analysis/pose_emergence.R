##Load R packages
library(tidyverse)
library(ggplot2)

##Import data
source("data_compiling/compile_demography.R")

#calculate standard error
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}
