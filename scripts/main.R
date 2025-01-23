
# always run before running any other files

# relevant packages
R_packages <- c("ggplot2", "rstan", "cowplot", "dplyr", "magrittr", "reshape2", "rockchalk", "here", "MASS",
                "microbiome", "phyloseq", "spatialEco", "loo", "splus2R", "diptest", "RColorBrewer", 
                "hilbertSimilarity", "modelbased", "ggpubr", "grid", "gridExtra", "ggeasy", "ggdist", "scales",
                "microViz", "parallel", "gganimate", "magick", "ggExtra", "purrr", "ggtext")

# download any new packages
new.packages <- R_packages[!(R_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load packages
for (p in R_packages) library(p, character.only = TRUE)

# get necessary functions
source(here("scripts", "functions", "functions.R"))

# compile model
EM_model <- stan_model(file = here("scripts", "model", "EM_model.stan"))
