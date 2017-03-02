setwd("DownstreamAnalysisOfNetworkStructure")
#The original data for running this comes from "makeJSOncombinedMethod.R"
#you can load this data from here:
load("testRobustness.RData")
load("clustersForNetworkDrawing.RData")
#this contains:
#edges_networkSizes[[2]]: all edges in 5000 edge network
#allRNA: gene expression only for these genes


#set seed
set.seed(123)

#select 10 random sets, each containing 5 TFs (might eventually increase to 20, but it might take a really long time!)
tfs=as.character(unique(edges_networkSizes[[2]][,1]))
tfSamples=lapply(c(1:20), function(i){sample(tfs, 5)})

#select 5 random subsets of test set vs. training set (50 random time points?)
rnaseq_testSet=lapply(c(1:5), function(i){sample(c(1:length(rownames(allRNA))), 50)})
rnaseq_trainSet=lapply(c(1:5), function(i){which(! (c(1:length(rownames(allRNA))) %in% rnaseq_testSet[[i]]))})


library("randomForest")
#for each set of genes
rf_results=sapply(tfSamples, function(tfSet){ 
#for each target gene
     sapply(as.character(colnames(allRNA)), function(target){
#for each set of test vs. train data
            mean(sapply(c(1:5), function(i){
              a=randomForest(allRNA[rnaseq_trainSet[[i]],tfSet], y=allRNA[rnaseq_trainSet[[i]], target], xtest=allRNA[rnaseq_testSet[[i]],tfSet], ytest=allRNA[rnaseq_testSet[[i]], target])
#calculate mse on test data for rf
              cor(a$test$predicted,allRNA[rnaseq_testSet[[i]],target])
  
            }))
     }) 
  
})

save(rf_results, file="rf_results_2.RData")

rf_results_median=apply(rf_results, 1, function(i){median(i)})
rf_results_sd=apply(rf_results, 1, function(i){sd(i)})
rf_results_range=apply(rf_results, 1, function(i){max(i)-min(i)})
rf_results_max=apply(rf_results, 1, function(i){max(i)})
rf_results_min=apply(rf_results, 1, function(i){min(i)})



#order by number of edges leading to gene
numParents=sapply(names(rf_results_median), function(i){
  length(which(as.character(edges_networkSizes[[2]][,2])==i))
})

numChildren=sapply(names(rf_results_median), function(i){
  length(which(as.character(edges_networkSizes[[2]][,1])==i))
})

#boxplot(numParents, rf_results_median)


#Conclusion: genes with more 'parents' are the ones that are most redundantly regulated (but it is not related to the sd of teh correlation)
boxplot(lapply((0:5), function(j){rf_results_median[which(numParents==j)]}))
boxplot(lapply((0:5), function(j){rf_results_sd[which(numParents==j)]}))

pdf(paste("Figure_robustnessBoxPlots.pdf", sep=""), width=8.5, height=5,pointsize = 10) 
par(mfrow=c(2,2));
par(mar=c(4, 4, 3, 3)+0.1);
par(cex=0.9)


boxplot(lapply((0:5), function(j){rf_results_min[which(numParents==j)]}), xlab="# edges in predicted network", ylab="minimum R", col="grey", ylim=c(-0.1,1))
boxplot(lapply((0:5), function(j){rf_results_max[which(numParents==j)]}), xlab="# edges in predicted network", ylab="maximum R", col="grey", ylim=c(-0.1,1))

boxplot(
  lapply(c(5, 7, 4, 2, 6, 1, 3), function(j){rf_results_min[names(clustersForNetworkDrawing)[which(clustersForNetworkDrawing==j)]]}), ylim=c(-0.1,1), xlab="cluster", ylab="minimum R", col=c("lightgreen", "pink", "darkgreen", "cornflowerblue",  "red", "orange", "tan"))

boxplot(
  lapply(c(5, 7, 4, 2, 6, 1, 3), function(j){rf_results_max[names(clustersForNetworkDrawing)[which(clustersForNetworkDrawing==j)]]}), ylim=c(-0.1,1), xlab="cluster", ylab="maximum R", col=c("lightgreen", "pink", "darkgreen", "cornflowerblue", "red", "orange", "tan"))
dev.off()









