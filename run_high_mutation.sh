#!/bin/bash

LOCI=$(python -c 'print([1e-5]*20)' | tr -d ',[]')

/home/local/USDA-ARS/javier.tabima/Documents/zhian_simulations/simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 13 --murate 1e-3 $LOCI --nloc 21 --outfile high_mutation1 &
/home/local/USDA-ARS/javier.tabima/Documents/zhian_simulations/simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --murate 1e-3 $LOCI --nloc 21 --outfile high_mutation2 &
/home/local/USDA-ARS/javier.tabima/Documents/zhian_simulations/simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 13 --murate 1e-3 $LOCI --nloc 21 --outfile high_mutation3 &
/home/local/USDA-ARS/javier.tabima/Documents/zhian_simulations/simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --murate 1e-3 $LOCI --nloc 21 --outfile high_mutation4 &
/home/local/USDA-ARS/javier.tabima/Documents/zhian_simulations/simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 13 --murate 1e-3 $LOCI --nloc 21 --outfile high_mutation5 &
/home/local/USDA-ARS/javier.tabima/Documents/zhian_simulations/simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --murate 1e-3 $LOCI --nloc 21 --outfile high_mutation6 &
/home/local/USDA-ARS/javier.tabima/Documents/zhian_simulations/simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 13 --murate 1e-3 $LOCI --nloc 21 --outfile high_mutation7 &
/home/local/USDA-ARS/javier.tabima/Documents/zhian_simulations/simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --murate 1e-3 $LOCI --nloc 21 --outfile high_mutation8

/home/local/USDA-ARS/javier.tabima/Documents/zhian_simulations/organizing/dir2feather.py --prefix twenty_loci --group_by sex -z -o zip_high_mutation
