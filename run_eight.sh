#!/usr/bin/sh

python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt2 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt3 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt4 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt5 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt6 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt7 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt8
