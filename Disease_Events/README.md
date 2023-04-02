# Disease Events (Fig. 1)

System requirements:

-GCC (checked with Apple clang version 12.0.5 (clang-1205.0.22.9))

-R (checked with R version 3.6.3 (2022-05-25))

-R packages ggplot.

Once R and gcc are available, no explicit installation, or special hardware is needed. When running the master script as described below, C codes will be compiled as external libraries. R libraries will be installed automatically, although installation success depends upon R's version. If R packages are not installed this way, please install them manually.

Please note that:
-Folder's names should be preserved. Otherwise errors will be raised as scripts need the correct paths for working properly.

##########################################################################################

Executing the following command from this folder:

bash master_script.sh

will produce an example of a simulation of the control cohort in the clinical trial formally equivalent to the one published on the main text. This simulations accounts for the number and frequency of disease events in each country of the multi-centric trial.

##########################################################################################

The 1_cases_test.R script produce the RCT simulation and it load codes in C as external libraries.
Following the execution of 1_cases_test.R, 2_Comb_placebo_arm.R will be executed to produce the plots of the event distributions and the result for the complete placebo arm.

##########################################################################################

Expected output: as a final result of the simulation, the algorithm will write the output for each type of event and its frequency (in Outputs/ folder) along with the distribution of events in each country and combined in the complete placebo arm.

