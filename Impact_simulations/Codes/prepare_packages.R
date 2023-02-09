#!/usr/bin/Rscript

my_packages <- c("truncnorm", "ggplot2", "cowplot", "gridExtra", "minpack.lm", "nlsr", "iterators", "foreach", "doParallel", "dplyr")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages

if(length(not_installed)){
    print(c("Packages not installed but needed:", not_installed))
    print("Installing")
    install.packages(not_installed, repos = "http://cran.us.r-project.org")
}

