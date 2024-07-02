
# load required libraries
library(data.table)


args<-commandArgs(trailingOnly = TRUE) 
subset.id<-as.numeric(args[1])
sample.size<-as.numeric(args[2])
input.path<-args[3]

community.profile<-fread(input.path, data.table = FALSE)
samples<-colnames(community.profile[,-1])

print(paste("subsetting data to ", sample.size, "samples", sep = " "))

# subset data
samples.sub<-sample(samples, size = sample.size)  
community.profile_sub<-community.profile[,c(1,which(colnames(community.profile) %in% samples.sub))]

# make sure the first column is named "#OTU ID", according to FastSpar specifications
colnames(community.profile_sub)<-c("#OTU ID", colnames(community.profile_sub)[-1])

# remove non-present OTUs
remove<-(which(rowSums(community.profile_sub[,-1]) == 0))

if (length(remove) > 0){
  community.profile_sub<-community.profile_sub[-remove,]
}

# remove OTUs present in <= 20% of samples (low prevalence OTUs cause issues in FastSpar)
binary_community.profile_sub<-community.profile_sub[,-1]
binary_community.profile_sub[binary_community.profile_sub > 0]<-1
remove<-(which(rowSums(binary_community.profile_sub) <= (sample.size*0.2)))

if (length(remove) > 0){
  community.profile_sub<-community.profile_sub[-remove,]
}


write.table(community.profile_sub, paste("output/data_subsets/community.profile_sub_", subset.id, ".tsv", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)



