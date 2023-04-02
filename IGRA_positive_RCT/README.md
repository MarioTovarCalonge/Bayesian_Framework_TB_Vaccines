# IGRA_positive_analysis

System requirements:

-GCC (checked with Apple clang version 12.0.5 (clang-1205.0.22.9))

-R (checked with R version 3.6.3 (2022-05-25))

-R packages: "ggplot2", "cowplot", "fanplot", "gridExtra", "kdensity", "KernSmooth", "viridis".

Once R and gcc are available, no explicit installation, or special hardware is needed. When running the master script as described below, C codes will be compiled as external libraries. R libraries will be installed automatically, although installation success depends upon R's version. If R packages are not installed this way, please install them manually.

Please note that:
-Folder's names should be preserved. Otherwise errors will be raised as scripts need the correct paths for working properly.

##########################################################################################

Executing the following command from this folder:

bash master_script.sh

will produce an example of a simulation of a TB vaccine clinical trial formally equivalent to the one published on the main text. Also, it will infer the bayesian posterior of each vaccine model using the proposed bayesian framework.


##########################################################################################

The MLE_test.R produce the RCT simulation and it combines codes in C as external libraries.
Following the execution of 1_MLE_test.R, 2_cloud_plot.R, and 3_Ridge_plot_Likelihoods.R will be executed automatically and they will produce the cloud and probability distributions from which the results in the main text are derived, and the model-posterior ranking

##########################################################################################

Expected output: as a final result of the simulation+vaccine inference procedure, the algorithm will write the output for each vaccine model (in Outputs/ folder) along with the model posteriors and inferred epsilon values for each model in the whole simulation.

