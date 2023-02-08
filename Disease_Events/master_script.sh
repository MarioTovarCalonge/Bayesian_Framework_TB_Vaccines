#!/bin/bash

cd Codes/
Rscript prepare_packages.R
Rscript 1_cases_test.R
Rscript 2_Comb_placebo_arm.R
cd ../

