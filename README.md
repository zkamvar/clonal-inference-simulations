# Zhian's simuPOP simulation scripts

## Purpose

The purpose of these scripts is to simulate populations with varying levels of
sexual reproduction and eventually varying levels of different population
genetic strategies such as admixture. 


## Broad workflow

Due to the raw size of the data associated, this repository does not track the
data generated. The broad workflow was initially supposed to be straightforward,
but due to unforseen bugs/circumstances/opportunities, many analyses had to be
re-run, so caveats like two genomic data output directories due to consequences
of linkage exist. 

### Generation of populations and statistics

The data are simulated with simuPOP, transferred to feather files, and then
statistics are generated in parallel on the CGRB infrastructure, generating 
several binary \*.rda files for downstream analysis.

1. Simulate and save populations
    - Simulate Populations with simuPOP and save as \*.pop files with scripts in
      [`simulations/`](simulations/), usually facilitated with run_\*.sh scripts or run_\*.txt
      files and [`recovery_helper.sh`](recovery_helper.sh).
    - Run unfinished simulations to completion with [`simulations/recover_ssr.py`](simulations/recover_ssr.py)
      and [`recovery_helper.sh`](recovery_helper.sh)
    - Convert \*.pop files to feather format and zip with
      [`organizing/dir2feather.py`](organizing/dir2feather.py)
    - Synch zipped feather files to CGRB infrastructure
2. Sample populations and calculate statistics for microsatellite data
    - Read zipped feather file, subsample, and apply
      [`analysis/analyze_and_save_ia.R`](analysis/analyze_and_save_ia.R) with [`analysis/unzip_and_analyze.sh`](analysis/unzip_and_analyze.sh).
      Subsampled data are saved as \*.DATA.rda files and are then used
      for further analyses. **Output dir: rda\_files**
    - Gather genotypic diversity statistics on \*.DATA.rda files with
      [`analysis/analyze_diversity_table.R`](analysis/analyze_diversity_table.R) **Output dir: diversity\_rda\_files**
    - Gather allelic diversity statistics on \*.DATA.rda files with
      [`analysis/analyze_locus_table.R`](analysis/analyze_locus_table.R) **Output dir: locus\_rda\_files**
    - Gather locus contribution statistics on \*.DATA.rda files with
      [`analysis/analyze_locus_contribution.R`](analysis/analyze_locus_contribution.R) **Output dir: locus\_contribution\_rda_files**
3. Sample populations and calculate statistics for SNP (genomic) data
    - Read zipped feather file, subsample, and apply
      [`analysis/genomic_analyze_and_save_ia.R`](analysis/genomic_analyze_and_save_ia.R) with 
      [`analysis/unzip_and_analyze.sh`](analysis/unzip_and_analyze.sh), shuffling each locus independently and 
      additionally inserting missing data at 1%, 5% and 10% and assessing ia
      (no significance analysis). **Output dir: genomic\_rda\_files**
    - Re-analyze significance testing of ia with
      [`analysis/genomic_re_analyze_ia.R`](analysis/genomic_re_analyze_ia.R), shuffling linkage blocks with
      blocksize equal to 1000. **Output dir: genomic\_rda\_files\_redux**
4. Sync results folder to local workstation in the directory `../simulation_results`

### Exploratory data analysis

These analyses are facilitated via the *Rmd files in the [`analysis/`](analysis/) folder. 
There are data processing steps and analysis steps. All derivative data is saved
in `../simulation_results`. The data processing/cleaning steps will only work
with the specific data associated with this analysis. This is due to caveats 
noted in the scripts themselves.

1. Process SSR data with [`analysis/ssr_data_cleaning.Rmd`](analysis/ssr_data_cleaning.Rmd). [*Note: this uses a LOT of RAM*][ramtweet]
2. Process SNP data with [`analysis/genomic_data_processing.Rmd`](analysis/genomic_data_processing.Rmd).
3. Calculate ROC analysis with [`analysis/ROC_calculation.Rmd`](analysis/ROC_calculation.Rmd).
4. Analyze SNP (genomic) data with [`analysis/genomic_data_analysis.Rmd`](analysis/genomic_data_analysis.Rmd)
5. Analyze SSR data with [`analysis/post_analysis.Rmd`](analysis/post_analysis.Rmd)
6. Analyze results of ROC analysis with [`analysis/ROC_Curve.Rmd`](analysis/ROC_Curve.Rmd)


## Current software environment for simulations:

 - Python 3.4
 - simuPOP 1.1.7
 - numpy 1.11.0
 - feather-format 0.3.0

I used conda to install all of these. [Here's the intro to conda][conda]. It's
worth noting that I did have some trouble getting everything to work correctly
since simuPOP can't be installed with python 3.5 (for some weird reason, black
magic, maybe? I dunno). These are my steps:

1. Install conda (on my OSX, I used `brew install Caskroom/cask/anaconda`, but
   there are instructions on the site.)
2. Update my modules via conda, Downgrade python and install simuPOP:

```
conda update --all 
conda install python=3.4 
conda install -c https://conda.binstar.org/bpeng simuPOP
``` 

Finally, I made sure that numpy and feather-format were installed:

```
conda install numpy
pip install feather-format # for transfer between python and R
```

## R package

To facilitate analysis of the simulations, I've written an R package called
"zksimanalysis". It requires a version of R >= 3.2.0, due to its dependency on
dplyr.

It can be installed from the root directory in R with:

```r
devtools::install("zksimanalysis")
```

### Installing on the CGRB cluster (CentOS 6.6)

The default version of g++ on the CGRB is 4.4, which is not great since both
*feather* and *vcfR* require C++11, which is only available in g++ version >= 
4.7. Luckily, there is a workaround. It turns out that there are developer tools
installed: http://superuser.com/a/542091

Matthew Peterson, (whom I owe a :beer:) suggested this to utilize the correct versions:

```sh
$ bash; source /opt/centos/devtoolset-1.1/enable; gcc --version
# gcc (GCC) 4.7.2 20121015 (Red Hat 4.7.2-5)
# Copyright (C) 2012 Free Software Foundation, Inc.
# This is free software; see the source for copying conditions.  There is NO
# warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

So, in order to install this on the CGRB, use:


```sh
SGE_Batch -r install_zksimanalysis -c 'bash; source /opt/centos/devtoolset-1.1/enable; R -e "devtools::install(\"zhian_simulations/zksimanalysis\")"'
```

## Current workflow/implementation

(2016-08-18)

Currently, the script can be executed from the top-level directory. I still need
to set this up so it dumps the results in a "results" directory, but that will
come later as I need it. All of the simulation scripts will live in the 
[`simulations/`][simulations] directory.

With default values:

```sh
python simulations/simulate_sex_rate.py
```

With a configuration file called `CONFIG.args`:

```sh
python simulations/simulate_sex_rate.py @CONFIG.args
```

Once all of these are created, they can be gathered up into the feather format:

```sh
python organizing/dir2feather.py --regex gen_10 --prefix foo --group_by sex --zip --out pillow
```

## [Modules][modules]

The `modules/` directory contains python scripts that can be loaded
as modules for classes and functions.

#### [random_alleles.py][random_alleles]

This module implements two classes, 

 - *zk_locus* will initialize a locus with a mutation rate, allele frequencies,
   repeat length, and allele names.
 - *zk_loci* contains a whole bunch of *zk_locus* classes. This will provide an
   accessor to all of the loci names, allele names, allele frequencies, etc.

#### [zk_utils.py][zk_utils]

This contains misc functions for evolution


#### [pop2df.py][pop2df]

This contains the code for converting population objects to pandas data frames
with the function `pops2df()`.


## [Simulation Scripts][simulations]

The following scripts are used for the simulations. Again, they are to be run 
from the top-level directory.

#### simulate_sex_rate.py

This will simulate populations with varying levels of sexual reproduction.
Currently, there is no population structure for the populations, but this can
change.


## [Organization][organizing]

This directory contains scripts used for organizing and transferring the output.

#### [dir2feather.py][dir2feather]

This will trawl a directory and convert the populations to the [feather
format][feather], which allows crosstalk between python and R for faster
read/write than traditional csv.

[conda]: http://conda.pydata.org/docs/intro.html
[modules]: ./modules
[simulations]: ./simulations
[random_alleles]: ./modules/random_alleles.py
[zk_utils]: ./modules/zk_utils.py
[pop2df]: ./modules/pop2df.py
[dir2feather]: ./organizing/dir2feather.py
[organizing]: ./organizing
[feather]: https://blog.rstudio.org/2016/03/29/feather/
[ramtweet]: https://twitter.com/ZKamvar/status/787777914288812033