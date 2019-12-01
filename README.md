# gwc
Method for generating genome-wide constraint metrics comparing observed to expected variation. Expected is calculated using an updated mutation rate model with sequence context + 32 additional features.

# Getting Started

This works well on R 3.2.2. It may work on later or older versions as well, but has not been thoroughly tested.

There are a number of R libraries required including dplyr, plyr, data.table, but nothing too unexpected.

The key script in the project is gwc/R/selection2.phylop.R - the most important inputs to this script are:
- tab delimited file defining elements for which elements we want to generate observed and expected variants
- config file describing the whole genome sequencing input(s) that will be used to generate the observed and expected variants

There is an example config file in gwc/R/gnomad_BRIDGE_config.R from running this script on the BRIDGE and gnomAD whole genomes.

Most of the required data should be available within this repo to run ```gwc/scripts/wg_2k_constraint_only``` - this takes a set of genome wide 2kb windows which has been pre-annotated with mutation rates and coverage and generates observed and expected variants.

One piece of data that is _not_ present because it is very large is the TSV files for the gnomAD whole genomes. I can provide this if needed, just email me! It is also on the gnomAD ftp site.

In order to re-run this on a different data set, a set of models describing the expected number of variants in the absence of selection must be produced using synonymous variants from the whole genome sequence data of interest. ```gwc/scripts/fit_synonymous_models_gnomAD_slurm``` can be adapted to do this for another project (e.g. Genomics England 100k genomes project).

