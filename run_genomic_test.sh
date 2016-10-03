#!/usr/bin/env bash

python simulations/simulate_sex_rate_gbs.py --POPSIZE 10000 --nseed 1 --nloc 10000 --outfile genomic1 &
python simulations/simulate_sex_rate_gbs.py --POPSIZE 10000 --nseed 1 --nloc 10000 --outfile genomic2 &
python simulations/simulate_sex_rate_gbs.py --POPSIZE 10000 --nseed 1 --nloc 10000 --outfile genomic3 &
python simulations/simulate_sex_rate_gbs.py --POPSIZE 10000 --nseed 1 --nloc 10000 --outfile genomic4 &
python simulations/simulate_sex_rate_gbs.py --POPSIZE 10000 --nseed 1 --nloc 10000 --outfile genomic5 &
python simulations/simulate_sex_rate_gbs.py --POPSIZE 10000 --nseed 1 --nloc 10000 --outfile genomic6 &
python simulations/simulate_sex_rate_gbs.py --POPSIZE 10000 --nseed 1 --nloc 10000 --outfile genomic7 &
python simulations/simulate_sex_rate_gbs.py --POPSIZE 10000 --nseed 1 --nloc 10000 --outfile genomic8