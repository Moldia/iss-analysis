# Title:        proximity_analysis_main.R
#
# Description:  Main script containing different high level procedures for
#                 visualization
#
#               Contains examples of a nearest neighbor analysis and some 
#                 plotting and explicit examples.
#                 
# Usage:        Make sure that source("proximity_analysis_source.R") is done
#
# Version:      2014-07-17 16:40:47 CEST
# Author:       Olle Nordesjö (olle.nordesjo[at]gmail.com)


#### --- Perform proximity analysis --- ###
# Below is a typical example of a complete proximity analysis


#--- Set parameters  --------------
par.workingdirectory="E:/Analysis_utilities/R/Poximity"
setwd(par.workingdirectory)
par.inputfile="proximity_analysis_sample_data_A549_celline.csv"  # Example file included
#par.subset=c("ALDH1A","RNF43","RAC1","MUC2","CEACAM1","AXIN2","KRT20","EPCAM")
par.subset=c("NNNN","ACTB")
par.distance=10 #in pixels

# First, import the data and subset it. Subset with inlude=T to include or
#   include=F to exclude

data.raw<-read.csv(file=par.inputfile)                                          
#data.subset<-SubsetData(data.raw,subset=par.subset,include=T)
data.subset<-subsetData(data.raw,subset=par.subset,include=F) 

# Inspect the plot to see if the coordinates seem to be correct
data.subset_plotXY<-plotXY(data.subset,size=0.1,alpha=0.9,facet=T);   
data.subset_plotXY                                                       

# Create a hash to generate identities to the transcript names
hash.names<-hash(names(table(droplevels(data.subset$name))),
                 as.numeric(factor(names(table(droplevels(data.subset$name))))))

# Generate matrix with replicates to perform bootstrapping
uNames<-unique(data.subset$name)                                    
nNames<-length(uNames)                                              
proximityReplicates<-array(0,dim=c(nNames,nNames,100))

# Generate all the randomized replicates
for (times in 1:100){
  print(times)
  data.rand<-generateRandomData(data.subset)
  pairwiseData<-calculateProximity(data.rand,par.distance)
  
  # loop through the pairwise distances and save the score in the matrix
  for(iteration in seq_along(pairwiseData[,1])){
    q.id=hash.names[[as.character(pairwiseData$query[iteration])]]
    t.id=hash.names[[as.character(pairwiseData$target[iteration])]]
    proximityReplicates[q.id,t.id,times]<-pairwiseData[iteration,]$score
  }  
}


# Now load the "real" data
data.process<-data.subset
pairwiseData<-calculateProximity(data.process,par.distance)

# Reordering of the levels so they reflect the frequency of occurrence
# table automatically sorts alphabetically. The "order" function reorders after frequency instead
names.sorted<-names(rev(table(droplevels(data.process$name))[order(table(droplevels(data.process$name)))]))
pairwiseData$query<-with(pairwiseData,
                              factor(query,levels=names.sorted))
pairwiseData$target<-with(pairwiseData,
                               factor(target,levels=names.sorted))


# Calculate some statistics for the pairwise comparisons
#   probability:  the fraction of times that the real proximityscore is
#                 higher than a randomly generated one.
#   deviation:   The real proximity score divided by the 
#                 median proximity score. (*100)

for (i in 1:length(pairwiseData[,1])){
  print(i)
  pairwiseData$probability[i]=calculateProbability(proximityReplicates, # Skip this, probably
                                                        pairwiseData$query[i],
                                                        pairwiseData$target[i])
  pairwiseData$deviation[i]=calculateDeviation(proximityReplicates,
                                                      pairwiseData$query[i],
                                                      pairwiseData$target[i])
}

# Generate a plot over how often transcripts occur near each other
#   Parameters are adjustable:
#     colType:   The variable that should be used to color the data points
#     sizeType:  The variable that should control the size of the data points.
#                 Remove if all should be of same size
#     colours:   Adjust according to taste

proximityPlot<-plotPairwise(pairwiseData,colType="probability",
                                              textsize=10)
proximityPlot<-plotPairwise(pairwiseData,colType="deviation",
                                     colours=c("red","red","grey","green","green"),
                                     textsize=10,sizeType="deviation")
proximityPlot

# Save the plot in a file
ggsave(filename=paste(par.workingdirectory,"proximityPlot",par.distance,
                      par.inputfile,".png",sep="_"),plot=proximityPlot,width=10,height=9)

# Calculate distance scores between transcripts based on the proximity profiles
#   Determine how similar two different transcripts are based on their 
#   proximity patterns. If both are near the same transcripts, this may yield 
#   information about their function. Refer to the function for details
similarityMatrix<-calculateProximitySimilarity(pairwiseData=pairwiseData)

# The heatmap shows more similar transcripts as hotter.
#   Only very obvious patterns should be noted, as the
#   strategy has not been demonstrated to be very sensitive.
plotSimilarityHeatmap(similarityMatrix)
