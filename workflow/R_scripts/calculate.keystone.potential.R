
# load required libraries
library(data.table)
library(reshape2)
library(igraph)
library(dplyr)


# build function to calculate relative node degree, closeness centrality and betweenness centrality
# correlation network keystones tend to have high closeness centrality
# however, if closeness centrality shows a strong positive linear correlation with degree, we recommend to not consider it

# get number of iterations from Rscript input
args<-commandArgs(trailingOnly = TRUE)

iterations<-as.numeric(args[1])
print(paste(iterations, "iterations", sep = " "))

sig.cutoff<-as.numeric(args[2])
print(paste("significant cutoff:",sig.cutoff, sep = " "))

use.closeness.centrality<-args[3]

centrality.norm<-function(graph,type=c("degree","closeness","betweenness"),centralization=FALSE)
{
  result<-NA
  g<-graph
  cent<-centralization
  if (!is.igraph(g)) {stop("Not a graph object")}
  if (type[[1]] == "degree") {
    if (!cent) result <- degree(g,mode="all")/(vcount(g)-1)
    else result <- (sum(max(degree(g,mode="all"))-degree(g,mode="all")))/((vcount(g)-1)*(vcount(g)-2))}
  else if (type[[1]] == "betweenness") {
    temp <- 2*betweenness(g,directed=F)/((vcount(g)-1)*(vcount(g)-2))
    if (!cent) result <- temp
    else  result <- sum(max(temp)-temp)/(vcount(g)-1)}
  else if (type[[1]] == "closeness") {
    if (!cent) result <- closeness(g,mode="all")
    else result <- (2*vcount(g)-3)*(sum(max(closeness(g,mode="all"))-closeness(g,mode="all")))/((vcount(g)-1)*(vcount(g)-2))}
  else {stop("this type is unavailable or mispelled")}
  return(result)
}


# create empty data frame to collect all correlations
correlations<-data.frame(source=NA, target=NA, correlation.strength=NA)

# add iteration
for (i in 1:iterations){

  # read in correlation data frame
  corr<-fread(paste("output/median_correlations/median_correlation_",i, ".tsv", sep = ""), data.table = FALSE) 
  # read in corresponding p-values data frame
  pvalues<-fread(paste("output/pvalues/pvalues_", i, ".tsv", sep = ""), data.table = FALSE) 
  
  # remove first column of correlation data frame containing OTU identifiers and instead set row names of data frame to these OTU identifiers
  r.names<-corr$`#OTU ID` 
  corr<-corr[,-1]
  rownames(corr)<-r.names
  
  # set upper triangular to 0 (as we are dealing with undirected correlations, upper and lower triangular are identical)
  corr[upper.tri(corr)]<-0
  
  # set diagonal to 0 (remove self-loops from correlations)
  diag(corr)<-0
  
  # melt data frame into long form and rename columns
  corr<-reshape2::melt(as.matrix(corr)) 
  colnames(corr)<-c("source", "target", "correlation.strength")
  
  # repeat with pvalues data frame
  pvalues<-pvalues[,-1]
  rownames(pvalues)<-r.names
  pvalues<-reshape2::melt(as.matrix(pvalues))
  colnames(pvalues)<-c("source", "target", "pvalue")
  
  # identify non-significant p-values and set corresponding correlations to 0
  not.sig<-which(pvalues$pvalue > sig.cutoff)
  corr$correlation.strength[not.sig]<-0
  
  # remove all 0s from list of iterations
  corr<-corr[corr$correlation.strength != 0,]
  
  # add significant correlations to list of all correlations
  correlations<-rbind(correlations, corr)
}

# remove first row containing NAs
correlations<-correlations[-1,]

# combine source and target column to species.pair column
correlations$species.pair<-apply(correlations[,c(1,2)], 1, paste0, collapse = "_:_" )

# calculate mean correlation strength across all identical species pairs
correlations_means<-
  correlations %>%
  group_by(species.pair) %>%
  summarise(mean.correlation.strength = mean(correlation.strength, na.rm = TRUE)) 

# count the occurrence of unique species pairs (= number of significant iterations) in list of significant correlations
correlations_counts<-as.data.frame(table(correlations$species.pair))
colnames(correlations_counts)<-c("species.pair", "significant.iterations")  
correlations_counts$species.pairs<-as.vector(correlations_counts$species.pairs)

# put together a data frame with source, target, mean correlation strength and number of significant iterations
correlation_pairs<-unique(correlations[,c(1,2,4)])
correlation.network<-merge(correlation_pairs, correlations_means)
correlation.network<-merge(correlation.network, correlations_counts)
correlation.network<-correlation.network[,-1]

# set mean correlation strength to 0 if the number of significant iterations is < 20% of iterations
correlation.network$mean.correlation.strength[correlation.network$significant.iterations < 2]<-0

# pull out source and target names
nameVals <- sort(unique(unlist(correlation.network[1:2])))

# construct 0 matrix of correct dimensions with row and column names
corr_adj <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))

# fill in the matrix with matrix indexing on row and column names
corr_adj[as.matrix(correlation.network[c("source", "target")])] <- correlation.network[["mean.correlation.strength"]]

# make it an adjacency matrix (all non-zero correlations are set to 1)
corr_adj<-as.matrix(ifelse (corr_adj!=0 ,1,0))

# build an undirected graph from adjacency matrix
g_total<-graph.adjacency(corr_adj, mode = "undirected")

# remove isolated vertices from the graph
isolated <- which(degree(g_total)==0)
g <- delete.vertices(g_total, isolated)

# calculate graph (=network) properties
degree<-centrality.norm(g,type="degree",centralization=F)
cb<-centrality.norm(g,type="betweenness",centralization=F)
cc<-centrality.norm(g,type="closeness",centralization=F)
tr<-transitivity(g,type="local")

# save network properties of each OTU in a data frame
stats<-data.frame(OTU=names(degree), Degree=degree, betweenness.centrality=cb, closeness.centrality=cc, transitivity=tr)
rownames(stats)<-NULL

# calculate keystone potential per OTU ((Degree * Transitivity) / Betweenness centrality)
# depending on the relationship between Degree and Closeness centrality (see comment above), you may change to:
# ((Degree * Transitivity * Closeness centrality) / Betweenness centrality)

if (use.closeness.centrality == "NO"){
  stats$keystone.potential<-((as.numeric(stats$Degree) *
                                as.numeric(stats$transitivity))/
                               as.numeric(stats$betweenness.centrality))
}

if (use.closeness.centrality == "YES"){
  stats$keystone.potential<-((as.numeric(stats$Degree) *
                                as.numeric(stats$transitivity) *
                                as.numeric(stats$closeness.centrality))/
                               as.numeric(stats$betweenness.centrality))
}

# betweenness centrality and transitivity may be 0, resulting in Inf or NaN values
# remove them from the output
stats$keystone.potential[stats$keystone.potential == "Inf"]<-0
stats$keystone.potential[stats$keystone.potential == "NaN"]<-0

write.table(stats, "output/keystone.potential.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



