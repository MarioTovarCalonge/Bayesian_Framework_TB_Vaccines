#!/bin/bash

cd Codes/
Rscript prepare_packages.R
Rscript 1_MLE_test.R
Rscript 2_Cloud_plot.R
cd ../

