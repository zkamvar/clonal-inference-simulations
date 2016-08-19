#!/usr/bin/sh

python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --output pt1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --output pt2 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --output pt3 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --output pt4 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --output pt5 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --output pt6 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --output pt7 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 12 --output pt8
