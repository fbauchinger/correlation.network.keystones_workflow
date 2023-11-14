# correlation.network.keystones_workflow

published in: TBD  
workflow for computing significant and robust correlation networks from subsampled community profiles with FastSpar 
  
for more information on FastSpar see:  
Watts SC, Ritchie SC, Inouye M, Holt KE. 2019. FastSpar: rapid and scalable correlation estimation for compositional data. Bioinformatics 35:1064–1066.
https://doi.org/10.1093/bioinformatics/bty734  
Github: https://github.com/scwatts/fastspar  
  
workflow.R  
R Script containing basic FastSpar usage and workflow for computing a robust consensus network from multiple FastSpar runs  
  
test_data:  
community_profile.tsv ... example input data for FastSpar (containing total reads; OTUs are rows; samples are columns)  
pvalues_[1:10].tsv ... example output files containing p-values of 10 different FastSpar runs (each run with a random subsample from a bigger dataset)  
observed_correlations_[1:10].tsv ... example output files containing observed correlations of 10 different FastSpar runs (each run with a random subsample from a bigger dataset)
