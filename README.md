# correlation.network.keystones_workflow

published in: TBD  
workflow for computing significant and robust correlation networks from subsampled community profiles with FastSpar 
  
for more information on FastSpar see:  
Watts SC, Ritchie SC, Inouye M, Holt KE. 2019. FastSpar: rapid and scalable correlation estimation for compositional data. Bioinformatics 35:1064–1066.
https://doi.org/10.1093/bioinformatics/bty734  
Github: https://github.com/scwatts/fastspar  

FastSpar can be installed through conda: 
conda install -c bioconda -c conda-forge fastspar
for installation from source see https://github.com/scwatts/fastspar

required R packages:
  - data.table
  - reshape2
  - igraph
  - dplyr

The folder "workflow" contains a test community profile (community_profile.tsv in "input" directory), R scripts to subset the community profile (subset_data.R) and to calculate keystone potential (calculate_keystone.potential.R) in "R_scripts" directory and a wrapper script (keystone.potential.wrapper.sh) combining both scripts with FastSpar
  
Usage:
  - download the entire "workflow" folder
  - run "bash keystone.potential.wrapper.sh" from within the "workflow" folder in a terminal to run the test case

Notes:
  - variables for FastSpar and subsampling can be set within the wrapper script (significance cutoff should be decreased and iterations should be increased for a real use case - see variables indicated with "TEST CASE")
  - make sure the community profile is formated as follows: OTUs are rows, samples are columns, row containing the OTU tags is named "#OTU ID"
  - by default, closeness centrality is not used to calculate keystone potential, as it tends to strongly positively correlate with degree
    -> this can be changed within the wrapper script, but will decrease keystone potential by a factor of ~ 100 (closeness centrality values are < 1)
    
log and error messages are provided in output/log/keystone.potential.wrapper_logs.txt


Troubleshooting:
  - if libmkl_rt.so.2 can't be found: conda install mkl
  - if "parallel fastspar" does not work: apt install parallel
  - "fastspar" command not found: activate conda environment in your terminal  
