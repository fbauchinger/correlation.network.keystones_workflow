#!/bin/bash

########################################################
#
# Wrapper for FastSpar and keystone potential analysis
#
# Usage:
# bash keystone.potential.wrapper.sh
#
# Author: Franziska Bauchinger
#
########################################################

########################################################
#
# FastSpar can be installed through conda or from source:
# https://github.com/scwatts/fastspar
#
########################################################


# add in your input path 
INPUT_PATH="input/community_profile.tsv" # OTUs are rows, samples are columns

# add your output path (directories will be created below)
OUTPUT_PATH="output"

# set subsampling variables
SUBSAMPLE_SIZE=50
ITERATIONS=10 # TEST CASE -> should be at least 1000

# set FastSpar variables
FASTSPAR_ITERATIONS=10 # TEST CASE -> should be at least 100
EXCLUDE_ITERATIONS=20
THRESHOLD=0.1
BOOTSTRAPS=10 # TEST CASE -> should be at least 1000
PERMUTATIONS=10 # TEST CASE -> should be at least 1000

# set variables for keystone potential calculation
SIG_CUTOFF=0.1 # p-value cutoff for observed correlations  TEST CASE -> should be 0.05 or lower
USE_CLOSENESS_CENTRALITY=NO # set to "YES" if you want to include closeness centrality 

# make output and log folders (if they don't exist yet)
mkdir -p $OUTPUT_PATH
mkdir -p $OUTPUT_PATH/bootstrap_counts
mkdir -p $OUTPUT_PATH/bootstrap_correlations
mkdir -p $OUTPUT_PATH/data_subsets
mkdir -p $OUTPUT_PATH/pvalues
mkdir -p $OUTPUT_PATH/median_correlations
mkdir -p $OUTPUT_PATH/log

# specify log file
LOGFILE=$OUTPUT_PATH/log/keystone.potential.wrapper_logs.txt




for ((ID=1;ID<=ITERATIONS;ID++)); do

	echo "Iteration" $ID
	
	echo "subsampling data"

	# make random subsamples
	Rscript R_scripts/subset_data.R ${ID} $SUBSAMPLE_SIZE $INPUT_PATH >> "${LOGFILE}" 2>&1

	echo "running FastSpar"	
	
	# run fastspar with variables set above
	fastspar --iterations $FASTSPAR_ITERATIONS --exclude_iterations $EXCLUDE_ITERATIONS --threshold $THRESHOLD --otu_table $OUTPUT_PATH/data_subsets/community.profile_sub_${ID}.tsv --correlation $OUTPUT_PATH/median_correlations/median_correlation_${ID}.tsv --covariance $OUTPUT_PATH/median_correlations/median_covariance_${ID}.tsv >> "${LOGFILE}" 2>&1

	fastspar_bootstrap --otu_table $OUTPUT_PATH/data_subsets/community.profile_sub_${ID}.tsv --number $BOOTSTRAPS --prefix $OUTPUT_PATH/bootstrap_counts/bootstrap_iter_${ID} >> "${LOGFILE}" 2>&1

	parallel fastspar --iterations $FASTSPAR_ITERATIONS --exclude_iterations $EXCLUDE_ITERATIONS --threshold $THRESHOLD --otu_table {} --correlation $OUTPUT_PATH/bootstrap_correlations/cor_{/} --covariance $OUTPUT_PATH/bootstrap_correlations/cov_{/} -i 2 ::: $OUTPUT_PATH/bootstrap_counts/* >> "${LOGFILE}" 2>&1

	fastspar_pvalues --otu_table $OUTPUT_PATH/data_subsets/community.profile_sub_${ID}.tsv --correlation $OUTPUT_PATH/median_correlations/median_correlation_${ID}.tsv --prefix $OUTPUT_PATH/bootstrap_correlations/cor_bootstrap_iter_${ID} --permutations $PERMUTATIONS --outfile $OUTPUT_PATH/pvalues/pvalues_${ID}.tsv >> "${LOGFILE}" 2>&1


done



echo "calculating keystone potential"

Rscript R_scripts/calculate_keystone.potential.R $ITERATIONS $SIG_CUTOFF $USE_CLOSENESS_CENTRALITY >> "${LOGFILE}" 2>&1

echo "Done"

echo "log and error file: output/log/keystone.potential.wrapper_logs.txt"




########################################################
#
# Troubleshooting
# if libmkl_rt.so.2 can't be found: conda install mkl
#
# if "paralle fastspar" does not work: apt install parallel
#
# "fastspar" command not found: activate conda environment in the terminal
#
########################################################






