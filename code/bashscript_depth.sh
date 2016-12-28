#!/bin/bash -e

R CMD BATCH '--args 5 20 25 1 99 50' rscript_depth.R

# Parameters
# 1 field number
# 2 lower limit of depth
# 3 upper limit of depth
# 4 relative depth
# 5 number of iterations
# 6 upper integration limit for calculation of p values