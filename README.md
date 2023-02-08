# Bayesian_Framework_TB_Vaccines

System requirements:

-GCC (checked with Apple clang version 12.0.5 (clang-1205.0.22.9))

-R (checked with R version 3.6.3 (2023-02-08))

-R packages: "ggplot2", "cowplot", "fanplot", "gridExtra", "kdensity", "KernSmooth", "viridis", "truncnorm", "minpack.lm", "nlsr", "iterators", "foreach", "doParallel", "dplyr".

Once R and gcc are available, no explicit installation, or special hardware is needed. When running the master scripts provided, C codes will be compiled as external libraries. R libraries will be installed automatically, although installation success depends upon R's version. If R packages are not installed this way, please install them manually.

Please note that:
-Folder's names should be preserved. Otherwise errors will be raised as scripts need the correct paths for working properly.

The code is structured in three directories, each one with a master script that executes all relevant codes and produces formally equivalent outcomes to those in the main text.

##########################################################################################
# Disease Events

Executing the master script produces an example of a simulation of the control cohort in the clinical trial. This simulations accounts for the number and frequency of disease events in each country of the multi-centric trial.

##########################################################################################
# IGRA_positive_analysis

Executing the master script produces an example of a simulation of a TB vaccine clinical trial. Also, it will infer the bayesian posterior of each vaccine model using the proposed bayesian framework.

##########################################################################################
# Impact simulations

Executing the master script produces an example of an impact forecast simulation formally equivalent to the ones published on the main text. This impact forecast represents the impact of each one of the seven proposed vaccine models along with the bayesian average that is proposed withing this framework.

