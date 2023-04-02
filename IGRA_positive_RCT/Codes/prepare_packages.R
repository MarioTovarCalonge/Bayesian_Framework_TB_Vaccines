#!/usr/bin/Rscript

my_packages <- c("ggplot2", "cowplot", "fanplot", "gridExtra", "kdensity", "KernSmooth", "viridis")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages

if(length(not_installed)){
    print(c("Packages not installed but needed:", not_installed))
    print("Installing")
    install.packages(not_installed, repos = "http://cran.us.r-project.org")
}

#Create Folder structure.
output_main_dir <- c("Outputs/")
if (!dir.exists(output_main_dir)){
dir.create(output_main_dir)
} else {
    print(paste0(output_main_dir, " already exists."))
}

for(m in 1:7){
    output_dir <- paste0(output_main_dir, "Model_", m)
    if (!dir.exists(output_dir)){
    dir.create(output_dir)
    } else {
        print(paste0(output_dir, " already exists."))
    }
}
