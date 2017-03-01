#!/bin/bash

# set which parameters to look at
declare -a params
params=("depth" "rel_force_mean")

function overlay() {
    # runs overlay_pvals.r on pvalue files the parameter specified
    declare -a files  # initialise 'files' as an array
    for i in {1..5}  # assumes 5 fields
    do
	# get start of file name
	# NOTE: latter part of file name can vary depending on parameters
	# passed to original program
	name_start="pvalues_"$1"_"$i
	# get only the first result from ls in case there's more than one
	files[i]=$(ls $name_start* | head -1)
    done
    # call overlay.pvals for the parameter specified
    R CMD BATCH "--args ${files[*]} $1" overlay_pvals.r
}

for param in ${params[*]}
do
    overlay $param
    
done
