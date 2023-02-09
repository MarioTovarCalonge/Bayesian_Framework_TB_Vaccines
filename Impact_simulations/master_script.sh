#!/bin/bash

cd Codes/
Rscript prepare_packages.R
bash Launch500.sh Indonesia
bash Launch500.sh India
bash Launch500.sh Ethiopia
cd ../

cd Impact_Outputs/Plot_Codes/
Rscript Bayesian_Impact.R
Rscript Merge.R
Rscript Combined_Impact_plots.R
cd ../../
