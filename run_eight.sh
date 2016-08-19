#!/usr/bin/sh

python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt1 > pt1.log 2>&1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt2 > pt1.log 2>&1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt3 > pt1.log 2>&1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt4 > pt1.log 2>&1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt5 > pt1.log 2>&1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt6 > pt1.log 2>&1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt7 > pt1.log 2>&1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --outfile pt8 > pt1.log 2>&1
