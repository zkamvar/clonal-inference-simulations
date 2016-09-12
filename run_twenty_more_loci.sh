#!/usr/bin env sh

LOCI=$(python -c 'print([1e-5]*20)' | tr -d ',[]')

python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 9 --murate $LOCI --nloc 20 --outfile twenty_loci9 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 10 --murate $LOCI --nloc 20 --outfile twenty_loci10 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 9 --murate $LOCI --nloc 20 --outfile twenty_loci11 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 10 --murate $LOCI --nloc 20 --outfile twenty_loci12 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 9 --murate $LOCI --nloc 20 --outfile twenty_loci13 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 10 --murate $LOCI --nloc 20 --outfile twenty_loci14 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 9 --murate $LOCI --nloc 20 --outfile twenty_loci15 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 10 --murate $LOCI --nloc 20 --outfile twenty_loci16
