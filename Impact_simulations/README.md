# Impact simulations

System requirements:

-GCC (checked with Apple clang version 12.0.5 (clang-1205.0.22.9))

-R (checked with R version 3.6.3 (2022-05-25))

-R packages: "truncnorm", "ggplot2", "cowplot", "gridExtra", "minpack.lm", "nlsr", "iterators", "foreach", "doParallel", "dplyr".

Once R and gcc are available, no explicit installation, or special hardware is needed. When running the master script as described below, C codes will be compiled as external libraries. R libraries will be installed automatically, although installation success depends upon R's version. If R packages are not installed this way, please install them manually.

Please note that:
-Folder's names should be preserved. Otherwise errors will be raised as scripts need the correct paths for working properly.
-Uncompress all the .zip files included, as they contain all the needed inputs for the code to work along with the required folder structure. Errors will be raised otherwise.

##########################################################################################

Executing the following command from this folder:

bash master_script.sh

will produce an example of an impact forecast simulation formally equivalent to the ones published on the main text. This impact forecast represents the impact of each one of the seven proposed vaccine models along with the bayesian average that is proposed withing this framework.

##########################################################################################

The master script automatically calls and launch all required codes to produce the simulation. 
This is a complex code which is divided in several smaller scripts that work together to produce the simulation. In a nutshell, inside the Codes/ folder, the Launch500.sh bash script generates all the inputs for each vaccine model based on the in-silico RCT outcomes. Then, it executes the master500.R script which executes the spreading model and computes the impact. 

Following the previous steps, the impact outcomes are curated and combined using the bayesian posterior of each model within this framework, and the combined bayesian average is computed. Finally, the graphical outcome is produced.

##########################################################################################

Expected output: as a final result of the simulation, the algorithm will write the impact forecasted for each vaccine model (in Impact_Outputs/ folder) along with the bayesian average in three high-burden countries as a graphical outcome (All_Eff_Combined.pdf).

