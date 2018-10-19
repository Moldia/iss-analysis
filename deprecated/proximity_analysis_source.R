# Title:        proximity_analysis_source.R
#
# Description:  Contains all source code for the functions 
#                 used in the proximity analysis
#
# Usage:        run source("proximity_analysis_source.R") once or 
#                 press ctrl+shift+enter to load all the functions
#                 After this, open proximity_analysis_main.R for further 
#                 instructions
#
# Version:      2014-07-17 16:40:47 CEST
# Author:       Olle Nordesjö (olle.nordesjo[at]gmail.com)

#--- Set up the required packages  --------------

require(ggplot2)  # Package for creating the plots
require(RANN)     # Nearest neighbor analysis
require(reshape)  # Converting data between different formats
require(reshape2) # Converting data between different formats
require(RColorBrewer) # Setting color spaces
require(xtable)   # Generating tables
require(hash)     # Create hashes (pairs of keys and values)

subsetData<-function(data.orig,subset=c("NNNN"),include=F){
  #   Subsets a data set.
  #   This is for excluding unwanted transcripts from the downstream analysis
  #   
  #   Args: 
  #     data.orig:  The original, unsubsetted dataframe
  #     subset:     For example c("VGEF","RAC1") will subset on VGEF and RAC1.
  #     include:    If set to TRUE (default), will include the subset.
  #                 Otherwise, subset will be excluded              
  #
  #   Returns: 
  #     A subsetted data frame for downstream analysis
  
  indices.subset<-which(data.orig$name%in% subset)
  
  if(include==T){data.include<-data.orig[indices.subset,]
  } else {
    data.include<-data.orig[-indices.subset,]
  }  
  data.include<-data.include[which(data.include$name!="NNNN"),]
  return(data.include)
}
generateRandomData <- function (data=data.subset) {
  # Generate sample data for testing:
  # This function will randomize the names of the data
  #   which will act as a form of bootstrapping.
  data.rand<-data
  data.rand$name<-sample(data.rand$name)
  return(data.rand)
}
calculateProbability<-function(proximityReplicates,quer,targ){
  # Checks the fraction of times that the replicates (with generateRandomData)
  # generate a lower proximity score than the pairwiseData.
  #
  # Example: quer="RAC1", targ= "VEGF" will render the probability that 
  #   RAC1 and VEGF occurs near each other by random
  # 
  # Args:
  #   proximityReplicates:  Replicates of datasets with shuffled names
  #   quer: The query transcript name (e.g. "RAC1")
  #   targ: The target transcript name (e.g. "VEGF")
  #
  # Returns:
  #   The fraction of times that a randomized dataset will have lower proximity
  #   score
  
  q.id=hash.names[[as.character(quer)]]
  t.id=hash.names[[as.character(targ)]]
  
  if(length(na.omit(proximityReplicates[q.id,t.id,]))>2){
    ecdffunction<-ecdf(na.omit(proximityReplicates[q.id,t.id,]))
    return(ecdffunction(pairwiseData[which(pairwiseData$target==targ&pairwiseData$query==quer),]$score)
           )
  } else {
    return(0)
  }
}
calculateDeviation<-function(proximityReplicates,quer,targ){
  # Calculates the factor by which the transcripts occur with
  # compared to shuffled data
  #
  # Example: quer="RAC1", targ= "VEGF" will determine how common it is that a 
  #   randomized dataset will have a smaller proximity measure. 
  # 
  # Args:
  #   quer: The query transcript name (e.g. "RAC1")
  #   targ: The target transcript name (e.g. "VEGF")
  #
  # Returns:
  #   The probability that a randomized dataset will have less probability of proximity
  
  q.id=hash.names[[as.character(quer)]]
  t.id=hash.names[[as.character(targ)]]
  
  if(length(proximityReplicates[q.id,t.id,])>2){
    score=pairwiseData[which(pairwiseData$target==targ&pairwiseData$query==quer),]$score
    standardDev=sd(proximityReplicates[q.id,t.id,])
    mu=mean(proximityReplicates[q.id,t.id,])
    return(
      (score-mu)/standardDev
    )
  } else {
    return(0)
  }
  # TODO: figure out a better statistic (maybe)
}
calculateProximity<-function(data,distance){
  #   Calculates scores for pairwise transcripts proximity.
  #   Uses a nearest neighbor paradigm to register the  transcripts nearer than *distance* pixels to another.
  #   Assumes a mean cell-radius of *distance* pixels (ca *distance*/3 µm)
  #
  #   Utilizes the RANN-package (nearest neighbor tools)
  #   
  #   Args: 
  #     data: a dataframe from SubsetData
  #     distance: the maximum distance (in pixels) to accept as "proximity"
  #
  #   Returns: 
  #     A data frame with query transcript, target transcript, and the score
  
  
  hash.propensity<-hash(names(table(data$name)),as.numeric(table(data$name))/sum(as.numeric(table(data$name)))) # hash.propensity - transcript frequency
  namelist<-as.character(sort(unique(data$name)))
  proximity.data<-data.frame(query=c(),target=c(),score=c())
  
  # For each transcript name
  for (name in namelist){
    
    nnqueryids<-which(data$name==name)
    nndata<-data[-nnqueryids,]
    # Coordinates
    nnquery.xy<-data[nnqueryids,c("global_X_pos","global_Y_pos")]
    nndata.xy<-data[-nnqueryids,c("global_X_pos","global_Y_pos")]
    
    # Actual nearest neighbor calculation
    nnsd<-nn2(data = nndata.xy, query = nnquery.xy, 
              min(length(nndata.xy[,1]),15)) # Calculate the k nearest transcripts
    
    # Index-definitions
#     indices.fartherthan<-which(nnsd$nn.dists>distance)
#     nnsd$nn.dists[indices.fartherthan]<-NA
#     nnsd$nn.idx[indices.fartherthan]<-NA
    indices.closerthan<-which(nnsd$nn.dists<distance)
    
    if (is.factor(nndata$name)){
      close<-droplevels(nndata$name[nnsd$nn.idx[indices.closerthan]])
    } else {
      close<-nndata$name[nnsd$nn.idx[indices.closerthan]]
    }
    
    scores<-table(close) #/      (values(hash.propensity)[names(table(close))]*values(hash.propensity)[name])
    proximity.data<-rbind(proximity.data,data.frame(query=rep(name,length(scores)),
                                                    target=names(scores),
                                                    score=as.numeric(scores)))
  }

  
  proximity.data
    # TODO: add write to file!
}
calculateProximitySimilarity<-function(pairwiseData){
  #   Calculates similarity scores between transcripts and exports a heatmap
  #   Uses pairwise scoring of the transcripts using the probability
  #     score from pairwiseData. The output is a clustered heatmap showing
  #     the pairwise similarity score. 
  #   
  #   Scoring Scheme:
  #         If both are 1                    +1 
  #         If one is 1 and the other 0,     -1
  #         Otherwise (NA, etc.)              0
  #   
  #   Args: 
  #     pairwiseData: a dataframe containing probability from calculateProbability
  #
  #   Returns: 
  #     A similarity matrix (similarityMatrix) with pairwise similarity scores
  #       between transcripts. 
  #     Bigger value means that the transcript's proximity profiles are similar

  # Cast the probability values into a matrix
  probabilities<-round(acast(pairwiseData,target~query,value.var="probability"))
  dims<-dim(probabilities)[1]
  similarityScores<-matrix(0,nrow=dims,ncol=dims)
  
  # Calculate the the pairwise similarity for each gene
  for(iteration in 1:dims){
    for(jteration in iteration:dims){
      i=round(as.numeric(probabilities[iteration,]))
      j=round(as.numeric(probabilities[jteration,]))
      
      sameSites<-na.omit(i==1&j==1)
      diffSites<-na.omit((i==1&j==0)|(i==0&j==1))
      
      similarityScores[iteration,jteration]<-sum(abs(sameSites))-sum(abs(diffSites))
    } 
  }
  
  similarityScores.df<-as.data.frame(similarityScores)
  rownames(similarityScores.df)<-rownames(probabilities)
  colnames(similarityScores.df)<-rownames(probabilities)
  
  distanceMatrix<-max(similarityScores.df)-similarityScores.df
  
  distanceMatrix.sym<-matrix(-max(similarityScores.df),dims,dims)+t(distanceMatrix)+distanceMatrix
  similarityMatrix.sym<-max(distanceMatrix.sym,na.rm=T)-distanceMatrix.sym
  for (i in 1:dims){similarityMatrix.sym[i,i]=NA}
  return(similarityMatrix.sym)
}
plotSimilarityHeatmap<-function(distanceMatrix){
  return(heatmap(as.matrix(distanceMatrix),Colv="Rowv",symm=T,keep.dendro=T))
}
plotXY<-function(data,size=1,alpha=0.9,facet=T){
  #   Produces a faceted plot showing the locations of the transcripts
  #   
  #   Args: 
  #     data:   a dataframe from SubsetData
  #     size:   the size of the points
  #     alpha:  the opacity of the dots (0 to 1). Lower if data frame is dense
  #     facet:  whether to facet on transcript classes. Good for dense data
  #
  #   Returns: 
  #     A transcript-faceted plot of x and y position
  
  plotObject<-ggplot(data,aes(global_X_pos,global_Y_pos))+
    geom_point(size=size,alpha=alpha,aes(color=factor(name)))+
    guides(color=FALSE)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  if(facet==T){
    plotObject<-plotObject+facet_wrap(~name)
  }
  plotObject
  #ggsave(filename="facetPlot.png",plot=plotObject,width=12,height=12)  
}
plotPairwise<-function(pairwiseData,textsize=15,colType="score",sizeType=1,
                                colours=c("red","lightblue","lightblue","lightblue","black")){
  proximityPlot<-ggplot(pairwiseData,aes(query,target))+
    geom_point(aes_string(color=colType,size=sizeType))+
    scale_color_gradientn(colours=colours)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=textsize),
          axis.text.y = element_text(size=textsize))
  proximityPlot
  return(proximityPlot)
}
plotFrequency<-function(data){
  ggplot(data,aes(x=name))+geom_histogram()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
}