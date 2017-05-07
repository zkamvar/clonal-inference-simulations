Simulation Analysis Data
========================

This directory contains the results from the analysis of simulated population
data. It can be downloaded from the CGRB server with: 

    rsync --update -tavz -e "ssh -p XXXXX" --exclude 'twenty_loci_results' \
    kamvarz@files.cgrb.oregonstate.edu:\
    /nfs1/BPP/Grunwald_Lab/home/kamvarz/simulation_analysis/results/ .

Note that the `XXXXX` must be replaced with the proper port number.

The scripts in `analysis/` will read and write to this directory.

### Note

I started this project saving the results in a directory outside of this 
repository. I have since been creating soft-links to the files.
