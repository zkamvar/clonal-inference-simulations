#!/usr/bin env sh

LOCI=$(python -c 'print([1e-5]*20)' | tr -d ',[]')

python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 3 --murate $LOCI --nloc 20 --outfile twenty_loci1 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 3 --murate $LOCI --nloc 20 --outfile twenty_loci2 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 3 --murate $LOCI --nloc 20 --outfile twenty_loci3 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 3 --murate $LOCI --nloc 20 --outfile twenty_loci4 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 3 --murate $LOCI --nloc 20 --outfile twenty_loci5 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 3 --murate $LOCI --nloc 20 --outfile twenty_loci6 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 3 --murate $LOCI --nloc 20 --outfile twenty_loci7 &
python simulations/simulate_sex_rate.py --POPSIZE 10000 --nseed 3 --murate $LOCI --nloc 20 --outfile twenty_loci8
