# Bayesian_Framework_TB_Vaccines

System requirements:

-GCC (checked with Apple clang version 12.0.5 (clang-1205.0.22.9))

-R (checked with R version 3.6.3 (2022-05-25))

-R packages: fanplot, ggplot, cowplot, gridExtra, kDensity, KernSmooth, truncnorm, minpack.lm, nlsr, iterators, foreach, doParallel, dplyr.

Once R and gcc are available, no explicit installation, or special hardware is needed. When running the master script as described below, C codes will be compiled as external libraries. If R packages are not installed, please install them manually.

Please note that:
-Outputs should be uncompressed in order to provide the code the needed folder structure to save results.

##########################################################################################

Executing the following command from this folder:

Rscript MLE_test.R

will produce an example of a simulation of a TB vaccine clinical trial formally equivalent to the one published on the main text. Also, it will infer the likelihood of each vaccine model using the methodology explained there. 

PLEASE DON'T CHANGE ANY FOLDER NAME. Otherwise errors will be raised as scripts need the correct paths for properly working.

##########################################################################################

The MLE_test.R produce the RCT simulation and it combines codes in C as external libraries.
Following the execution of MLE_test.R, executing cloud_plot.R via Rscript will produce the cloud and probability distributions from which the results in the main text are derived.

##########################################################################################

Expected output: as a final result of the simulation+vaccine inference procedure, the algorithm will write the output for each vaccine model (in Outputs/ folder) along with the likelihood and inferred epsilon value for each model in the whole simulation.

