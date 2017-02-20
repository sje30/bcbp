#!/bin/bash -e

# Parameters
# 1 field number
# 2 number of iterations
# 3 upper integration limit for calculation of p values
# 4 parameter to be tested for sensitivity
# 5 lower limit of above parameter
# 6 upper limit of above parameter
# 7 increment between lower limit and upper limit

# parameter to be tested for sensitivity can be:
# 'depth', 'cone_excl', 'bc_excl_mean', 'bc_excl_sd', 'bc_excl_trunc', 'bcbp_excl_mean',
# 'bcbp_excl_sd', 'bcbp_excl_trunc', 'max_dendr_length', 'rel_force_mean'

# all combinations of field number and parameter to be tested are listed here with
# specific parameters. uncomment the line you want to run

# run this:
# -jN says to run N jobs in parallel.
# the sed command allows us to strip out the comments.
# sed '/^#/d' runs2.sh | nice -19 parallel -j40 


Rscript sensitivity_analysis.r 1 99 50 depth 20 25 1
Rscript sensitivity_analysis.r 2 99 50 depth 20 25 1
Rscript sensitivity_analysis.r 3 99 50 depth 20 25 1
Rscript sensitivity_analysis.r 4 99 50 depth 20 25 1
Rscript sensitivity_analysis.r 5 99 50 depth 20 25 1

Rscript sensitivity_analysis.r 1 99 50 cone_excl 10 13 0.5
Rscript sensitivity_analysis.r 2 99 50 cone_excl 10 13 0.5
Rscript sensitivity_analysis.r 3 99 50 cone_excl 10 13 0.5
Rscript sensitivity_analysis.r 4 99 50 cone_excl 10 13 0.5
Rscript sensitivity_analysis.r 5 99 50 cone_excl 10 13 0.5

Rscript sensitivity_analysis.r 1 99 50 bc_excl_mean 28 32 0.5
Rscript sensitivity_analysis.r 2 99 50 bc_excl_mean 36 39 0.5
Rscript sensitivity_analysis.r 3 99 50 bc_excl_mean 26 30 0.5
Rscript sensitivity_analysis.r 4 99 50 bc_excl_mean 31 35 0.5
Rscript sensitivity_analysis.r 5 99 50 bc_excl_mean 33 38 0.5

Rscript sensitivity_analysis.r 1 99 50 bc_excl_sd 4 8 0.5
Rscript sensitivity_analysis.r 2 99 50 bc_excl_sd 3 4.5 0.5
Rscript sensitivity_analysis.r 3 99 50 bc_excl_sd 2.5 4 0.5
Rscript sensitivity_analysis.r 4 99 50 bc_excl_sd 5 8 0.5
Rscript sensitivity_analysis.r 5 99 50 bc_excl_sd 5 8 0.5

Rscript sensitivity_analysis.r 1 99 50 bc_excl_trunc 15 19 0.5
Rscript sensitivity_analysis.r 2 99 50 bc_excl_trunc 30 33 0.5
Rscript sensitivity_analysis.r 3 99 50 bc_excl_trunc 21 23 0.5
Rscript sensitivity_analysis.r 4 99 50 bc_excl_trunc 16 20 0.5
Rscript sensitivity_analysis.r 5 99 50 bc_excl_trunc 21 23 0.5

Rscript sensitivity_analysis.r 1 99 50 bcbp_excl_mean 16 20 0.5
Rscript sensitivity_analysis.r 2 99 50 bcbp_excl_mean 20 24 0.5
Rscript sensitivity_analysis.r 3 99 50 bcbp_excl_mean 17 21 0.5
Rscript sensitivity_analysis.r 4 99 50 bcbp_excl_mean 16 20 0.5
Rscript sensitivity_analysis.r 5 99 50 bcbp_excl_mean 20 24 0.5

Rscript sensitivity_analysis.r 1 99 50 bcbp_excl_sd 3 7 0.5
Rscript sensitivity_analysis.r 2 99 50 bcbp_excl_sd 4 8 0.5
Rscript sensitivity_analysis.r 3 99 50 bcbp_excl_sd 3 7 0.5
Rscript sensitivity_analysis.r 4 99 50 bcbp_excl_sd 3 7 0.5
Rscript sensitivity_analysis.r 5 99 50 bcbp_excl_sd 6 10 0.5


Rscript sensitivity_analysis.r 1 99 50 bcbp_excl_trunc 7 9 0.5
Rscript sensitivity_analysis.r 2 99 50 bcbp_excl_trunc 6 8 0.5
Rscript sensitivity_analysis.r 3 99 50 bcbp_excl_trunc 4.5 6.5 0.5
Rscript sensitivity_analysis.r 4 99 50 bcbp_excl_trunc 4.5 6.5 0.5
Rscript sensitivity_analysis.r 5 99 50 bcbp_excl_trunc 4.5 6.5 0.5

Rscript sensitivity_analysis.r 1 99 50 max_dendr_length 41 47 0.5
Rscript sensitivity_analysis.r 2 99 50 max_dendr_length 41 47 0.5
Rscript sensitivity_analysis.r 3 99 50 max_dendr_length 41 47 0.5
Rscript sensitivity_analysis.r 4 99 50 max_dendr_length 41 47 0.5
Rscript sensitivity_analysis.r 5 99 50 max_dendr_length 41 47 0.5

Rscript sensitivity_analysis.r 1 99 50 rel_force_mean 0.65 0.85 0.05
Rscript sensitivity_analysis.r 2 99 50 rel_force_mean 0.80 1.00 0.05
Rscript sensitivity_analysis.r 3 99 50 rel_force_mean 0.70 0.90 0.05
Rscript sensitivity_analysis.r 4 99 50 rel_force_mean 0.75 0.85 0.05
Rscript sensitivity_analysis.r 5 99 50 rel_force_mean 0.85 1.05 0.05
