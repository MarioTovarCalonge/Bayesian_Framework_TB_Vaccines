#!/bin/bash

cd Codes/
Rscript prepare_packages.R
bash Launch500.sh Indonesia
bash Launch500.sh India
bash Launch500.sh Ethiopia
cd ../

cd Impact_Outputs/Plot_Codes/
Rscript 1_Bayesian_impact.R
Rscript 2_Merge.R
Rscript 3_Combined_Impact_plots.R
cd ../../
