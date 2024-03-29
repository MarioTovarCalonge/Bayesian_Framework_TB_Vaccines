#!/bin/bash

country=$1
aux1='M'

for i in `seq 1 7`;
do
    aux2=$(printf 'M%d' "$i")
    echo $aux2
    Rscript target_models_create_profile500.R $country $aux2 3
    Rscript master500.R $country 3
    Rscript Impact500_Target.R $country $i 3
done

for i in `seq 1 7`;
do
    aux2=$(printf 'M%d' "$i")
    echo $aux2
    Rscript target_models_create_profile500.R $country $aux2 0
    Rscript master500.R $country 3
    Rscript Impact500_Target.R $country $i 0
done
