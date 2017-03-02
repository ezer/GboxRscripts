#read in sequences
setwd("DAPseqAnalysis")
library("DNAshapeR")
library("vioplot")

seqCol0=read.table("geneTable")

####save it as fasta (for Homer2)
listForFasta=c()
for(i in c(1:length(seqCol0[,1]))){
  listForFasta=c(listForFasta, paste(">", as.character(seqCol0[i,1]), sep=""), as.character(seqCol0[i,6]))
}
write.table(listForFasta, file="GboxPromoters.fasta", row.names = FALSE, col.names = FALSE, quote=FALSE)

#find subset of genes with exactly one perfect G-box and extract flanking sequence
flank=sapply(c(1:length(seqCol0[,1])), function(i){
  spl=strsplit(as.character(seqCol0[i,6]), "CACGTG")
  if(length(spl[[1]])!=2){
    c(0,0)
  }else{
    
    c(substring(spl[[1]][1], nchar(spl[[1]][1])-9, nchar(spl[[1]][1])), 
      substring(spl[[1]][2], 1, 10))
  }
})
colnames(flank)=seqCol0[,1]
flank=flank[,which(flank[1,]!="0")]
flank=flank[,which(nchar(as.character(flank[1,]))==10 & nchar(as.character(flank[2,]))==10)]

#calculate parameters
sequenceParams=apply(flank, 2, function(i){
  c(unlist(lapply(c("A", "C", "G", "T"), function(base){
    a=lapply(c(1:10), function(j){
      if(substr(i[1], j, j)==base){1}else{0}
    })
    names(a)=paste(base, ":-", c(10:1), sep="")
    a
  })),
  
  unlist(lapply(c("A", "C", "G", "T"), function(base){
    a=lapply(c(1:10), function(j){
      if(substr(i[2], j, j)==base){1}else{0}
    })
    names(a)=paste(base, ":", c(1:10), sep="")
    a
  })))
  
  
  
})

input=apply(flank, 2, function(i){
  
  paste(i[1], "CACGTG", i[2], sep="")
  
})

listInput=c()
for(i in c(1:length(input))){
  listInput=c(listInput, paste(">", names(input)[i], sep=""), input[i])
}

write.table(listInput, row.names = FALSE, col.names = FALSE, quote=FALSE, file="input.txt")

shapeParams=getShape('input.txt')

shapeParams2=lapply(c(1:4), function(j){
  i=shapeParams[[j]]
  rownames(i)=names(input)
  colnames(i)=paste(names(shapeParams)[j], c(1:length(i[1,])), sep="")
  i[,which(apply(i, 2, function(t){!any(is.na(t))}))]
  
})

shapeParams3=cbind(shapeParams2[[1]], shapeParams2[[2]], shapeParams2[[3]], shapeParams2[[4]])
allParams=cbind(shapeParams3, t(sequenceParams))

#get clusters
cl=list()
cl[[1]]<-as.character(read.table('cl1.txt' ,header=FALSE, sep="\t")$V1)
cl[[2]]<-as.character(read.table('cl2.txt' ,header=FALSE, sep="\t")$V1)
cl[[3]]<-as.character(read.table('cl3.txt' ,header=FALSE, sep="\t")$V1)
cl[[4]]<-as.character(read.table('cl4.txt' ,header=FALSE, sep="\t")$V1)

allParams2=allParams[which(rownames(allParams) %in% c(cl[[1]], cl[[2]])),]

set.seed(123)
#test random forest: full model
accuracy=sapply(c(1:100), function(id){
  #training set
  trainIDs=sample(c(1:dim(allParams2)[1]), floor(0.66*dim(allParams2)[1]))
  trainSet=allParams2[trainIDs,]
  trainY=sapply(rownames(trainSet), function(i){ if(i %in% cl[[1]]){1}else{2}})
  #testing set
  testIDs=c(1:dim(allParams2)[1])[which(! c(1:dim(allParams2)[1]) %in% trainIDs)]
  testSet=allParams2[testIDs,]
  testY=sapply(rownames(testSet), function(i){ if(i %in% cl[[1]]){1}else{2}})
  
  #train rf
  #test accuracy in testing set
  
  length(which(randomForest(trainSet, as.factor(trainY), testSet, as.factor(testY))$test$predicted==testY))/length(testY)
})


#test random forest: shape only
accuracy_shape=sapply(c(1:100), function(id){
  #training set
  trainIDs=sample(c(1:dim(allParams2)[1]), floor(0.66*dim(allParams2)[1]))
  trainSet=allParams2[trainIDs,c(1:90)]
  trainY=sapply(rownames(trainSet), function(i){ if(i %in% cl[[1]]){1}else{2}})
  #testing set
  testIDs=c(1:dim(allParams2)[1])[which(! c(1:dim(allParams2)[1]) %in% trainIDs)]
  testSet=allParams2[testIDs,c(1:90)]
  testY=sapply(rownames(testSet), function(i){ if(i %in% cl[[1]]){1}else{2}})
  
  #train rf
  #test accuracy in testing set
  
  length(which(randomForest(trainSet, as.factor(trainY), testSet, as.factor(testY))$test$predicted==testY))/length(testY)
})


#test random forest: sequence only
accuracy_seq=sapply(c(1:100), function(id){
  #training set
  trainIDs=sample(c(1:dim(allParams2)[1]), floor(0.66*dim(allParams2)[1]))
  trainSet=allParams2[trainIDs,c(91:dim(allParams2)[2])]
  trainY=sapply(rownames(trainSet), function(i){ if(i %in% cl[[1]]){1}else{2}})
  #testing set
  testIDs=c(1:dim(allParams2)[1])[which(! c(1:dim(allParams2)[1]) %in% trainIDs)]
  testSet=allParams2[testIDs,c(91:dim(allParams2)[2])]
  testY=sapply(rownames(testSet), function(i){ if(i %in% cl[[1]]){1}else{2}})
  
  #train rf
  #test accuracy in testing set
  
  length(which(randomForest(trainSet, as.factor(trainY), testSet, as.factor(testY))$test$predicted==testY))/length(testY)
})



#random classifier
accuracy_random=sapply(c(1:100), function(id){
  #training set
  trainIDs=sample(c(1:dim(allParams2)[1]), floor(0.66*dim(allParams2)[1]))
  trainSet=allParams2[trainIDs,]
  trainY=sapply(rownames(trainSet), function(i){ if(i %in% cl[[1]]){1}else{2}})
  #testing set
  testIDs=c(1:dim(allParams2)[1])[which(! c(1:dim(allParams2)[1]) %in% trainIDs)]
  testSet=allParams2[testIDs,]
  testY=sapply(rownames(testSet), function(i){ if(i %in% cl[[1]]){1}else{2}})
  
  #train rf
  #test accuracy in testing set
  
  length(which(randomForest(trainSet, as.factor(trainY), 
                            testSet)$test$predicted==sample(testY, length(testY))))/length(testY)
})


vioplot(accuracy, accuracy_shape, accuracy_seq, accuracy_random, col="grey")

#best features of sequence model:
seqFeatures=importance(randomForest(allParams2[,c(91:dim(allParams2)[2])], as.factor(sapply(rownames(allParams2), function(i){ if(i %in% cl[[1]]){1}else{2}}))))
groups=sapply(c("A", "C", "G", "T"), function(i){
  which(substring(names(seqFeatures[,1]), 1,1)==i)})

seqFeaturesTable=cbind(seqFeatures[groups[,1]],
                       seqFeatures[groups[,2]],
                       seqFeatures[groups[,3]],
                       seqFeatures[groups[,4]])
rownames(seqFeaturesTable)=c(-10:-1, 1:10)
barplot(t(seqFeaturesTable), beside=TRUE, col=c("darkgreen", "darkblue", "goldenrod", "red"), ylab="importance (mean decrease Gini)", xlab="position relative to G-box")

#for weblogo: try to make motif of cl1 vs cl2
write.table(input[which(as.character(names(input)) %in% as.character(cl[[1]]))], file="cl1_weblogo", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(input[which(as.character(names(input)) %in% as.character(cl[[2]]))], file="cl2_weblogo", row.names = FALSE, col.names = FALSE, quote=FALSE)

#what percentage of cluster 1 have GA or TC?  
apply(allParams2[which(as.character(rownames(allParams2)) %in% as.character(cl[[1]])), c("A:-2", "C:-2", "G:-2", "T:-2", 
                                                                                         "A:-1", "C:-1", "G:-1", "T:-1",
                                                                                         "A:1", "C:1", "G:1", "T:1",
                                                                                         "A:2", "C:2", "G:2", "T:2")], 2, function(i){sum(i)})/length(cl[[1]])

apply(allParams2[which(as.character(rownames(allParams2)) %in% as.character(cl[[2]])), c("A:-2", "C:-2", "G:-2", "T:-2", 
                                                                                         "A:-1", "C:-1", "G:-1", "T:-1",
                                                                                         "A:1", "C:1", "G:1", "T:1",
                                                                                         "A:2", "C:2", "G:2", "T:2")], 2, function(i){sum(i)})/length(cl[[2]])



shapeFeatures=importance(randomForest(allParams2[,c(1:90)], as.factor(sapply(rownames(allParams2), function(i){ if(i %in% cl[[1]]){1}else{2}}))))
shapeFeaturesTable=cbind(c(shapeFeatures[c(1:22),], 0),
                         shapeFeatures[c(23:45),],
                         c(shapeFeatures[c(46:67),],0),
                         shapeFeatures[c(68:90),])
colnames(shapeFeaturesTable)=c("MGW", "HelT", "ProT", "Roll")


pdf(paste("Figure_randomForest.pdf", sep=""), width=4, height=4,pointsize = 10) 
par(mfrow=c(1,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=0.9)

vioplot(100*accuracy, 100*accuracy_seq, 100*accuracy_shape, 100*accuracy_random, col="grey", names=c("shape+seq", "seq", "shape", "random"), ylim=c(40, 100))

dev.off()



pdf(paste("Figure_randomForestVariableImportance.pdf", sep=""), width=7, height=5,pointsize = 10) 
par(mfrow=c(2,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=0.9)
barplot(t(seqFeaturesTable), beside=TRUE, col=c("darkgreen", "darkblue", "goldenrod", "red"), ylab="importance (mean decrease Gini)", xlab="position relative to G-box")
legend(1, 40, c("A", "C", "G", "T"), fill=c("darkgreen", "darkblue", "goldenrod", "red"))
barplot(shapeFeaturesTable, beside=TRUE, ylab="importance (mean degrease Gini)", xlab="DNA shape", col="grey")

dev.off()

