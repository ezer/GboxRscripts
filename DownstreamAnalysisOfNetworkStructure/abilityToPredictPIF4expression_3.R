#####Step 3: Random forest: can network predict pif4-101 gene expression better than a random network?

############for each gene: 
#find its #1 link
#find the number of links
setwd("../../DownstreamAnalysisOfNetworkStructure")
load("rna_seq_filtered.RData")
#Step4: Let's do this again, but for a random set of genes
#Step 1: get all the edges with a PIF4 (rank>10,000)
set.seed(302)
avg_rank_best=avg_rank[which(avg_rank<10000)]
edge_ends=sapply(names(avg_rank_best), function(i){substr(i, 11,19)})
edge_starts=sapply(names(avg_rank_best), function(i){substr(i, 1,9)})
unique_edge_ends=unique(edge_ends)
edges_with_pif4=edge_ends[which(as.character(edge_starts)==pif4)]


#step 2: for each of these, get a list of all the other genes that its dependent on
dependent_genes=lapply(unique_edge_ends, function(i){
  edge_starts[which(as.character(edge_ends)==as.character(i))]
})

#dependent_genes_all=dependent_genes
#dependent_genes=sample(dependent_genes,1000)

#Get how consistent the pattern in between Col-0 and pif4-101

unvariable=sapply(as.character(unique_edge_ends), function(i){
  pif4_id=which(as.character(rna_seq_filtered[["pif4_22"]][,1])==i)
  pif4_temp=c(rna_seq_filtered[["pif4_22"]][pif4_id,c(2,3,5:9)], rna_seq_filtered[["pif4_27"]][pif4_id,2:9])
  col0_id=which(as.character(rna_seq_filtered[["col0_22"]][,1])==i)
  col0_temp=c(rna_seq_filtered[["col0_22"]][col0_id, c(2,3,5:9)], rna_seq_filtered[["col0_27"]][col0_id,2:9])
  cor(as.numeric(col0_temp), as.numeric(pif4_temp))    
})


#rank_of_dependent_genes=avg_rank[names(dependent_genes)]
#dependent_genes=dependent_genes[order(rank_of_dependent_genes)]
#avg_rank[names(dependent_genes)]

#for each of these sets, train a rf model of gene expression
library("randomForest")
arr=read.table("gbox_minus_pif_data_set.txt")
colnames(arr)=as.character(t(arr[1,]))
arr=arr[-1,]
rf_models=lapply(dependent_genes, function(i){
  input=arr[,i]
  output=arr[,substr(names(i)[1], 11, 19)]
  randomForest(as.matrix(input), as.numeric(as.character(output)))
})

random_dependent_genes=lapply(dependent_genes, function(i){
  colnames(arr)[sample(c(1:length(colnames(arr))), length(i))]
})

rf_models_random=lapply(c(1:length(dependent_genes)), function(j){
  i=dependent_genes[[j]]
  input=arr[,random_dependent_genes[[j]]]
  output=arr[,substr(names(i)[1], 11, 19)]
  # print(output[1:4])
  randomForest(as.matrix(input), as.numeric(as.character(output)))
})


r_mat2=sapply(c(1:length(dependent_genes)), function(i){
  temp=substr(names(dependent_genes[[i]][1]),11,19)
  output=as.numeric(as.character(arr[,temp]))
  # print(output)
  #print(substr(names(dependent_genes[i])[1], 11, 19))
  c(cor(output, rf_models_random[[i]]$predicted),
    cor(output, rf_models[[i]]$predicted))
})


#Now, I can test the rf models on the pif4 ko data
#Now, I can test the rf models on the pif4 ko data
pif4_22=rna_seq_filtered[["pif4_22"]]
rownames(pif4_22)=pif4_22[,1]
pif4_27=rna_seq_filtered[["pif4_27"]]
rownames(pif4_27)=pif4_27[,1]

shared=rownames(pif4_22)[which(rownames(pif4_22) %in% rownames(pif4_27))]
all_pif4=t(cbind(pif4_22[shared,2:8], pif4_27[shared,2:9]))
all_pif4[,pif4]=rep(0, length(rownames(all_pif4)))
all_pif4=all_pif4[-3,]


pif4_rf_models=sapply(c(1:length(dependent_genes)), function(i){
  temp=dependent_genes[[i]]
  #temp=temp[which(temp!=pif4)]
  #print(temp)
  if(length(temp)!=1){
    input=all_pif4[,temp]
    #print(length(input))
    #print(dim(input))
    output=all_pif4[,substr(names(dependent_genes[[i]])[1], 11, 19)]
    p=predict(rf_models[[i]], input)
    cor(output, p)
  }else{
    input=as.matrix(all_pif4[,temp])
    #print(length(input))
    #print(dim(input))
    output=all_pif4[,substr(names(dependent_genes[[i]])[1], 11, 19)]
    p=predict(rf_models[[i]], input)
    cor(output, p)
  } })

pif4_rf_models_random=sapply(c(1:length(dependent_genes)), function(i){
  temp=random_dependent_genes[[i]]
  #temp=temp[which(temp!=pif4)]
  #print(temp)
  if(length(temp)!=1){
    input=all_pif4[,temp]
    #print(length(input))
    #print(dim(input))
    output=all_pif4[,substr(names(dependent_genes[[i]])[1], 11, 19)]
    p=predict(rf_models_random[[i]], input)
    cor(output, p)
  }else{
    input=as.matrix(all_pif4[,temp])
    #print(length(input))
    #print(dim(input))
    output=all_pif4[,substr(names(dependent_genes[[i]])[1], 11, 19)]
    p=predict(rf_models_random[[i]], input)
    cor(output, p)
  } })


r_mat_pif42=rbind(pif4_rf_models_random, pif4_rf_models)

pdf('Figure_abilityToPredictExpressionOfpif4_rnaseq_all.pdf', width=4.5, height=4.5 ,pointsize = 10)
par(mfrow=c(1,1));
par(mar=c(4, 4, 4, 4))
par(cex=1);
plot(r_mat_pif42[1,], r_mat_pif42[2,], xlab="random network", ylab="predicted network", pch=20, col=rgb(0.5, 0.5, 0.5, 0.15), xlim=c(-0.8, 1), ylim=c(-0.8, 1))
#points(r_mat_pif4[1,], r_mat_pif4[2,], col="red", pch=19)
#points(r_mat_pif4[1,], r_mat_pif4[3,], col="cornflowerblue", pch=19)
abline(c(0,1))
abline(v=0) 
abline(h=0)
#legend(0.1, -0.3, c("randomly selected genes with G-boxes", "genes with PIF4 as predicted target", "PIF4 removed from predicted network"), pch=c(20, 19, 19), col=c("grey", "red", "cornflowerblue"), cex=0.8, bty='n')

dev.off()

pdf('Figure_abilityToPredictExpressionOfpif4_rnaseq_all_CompareVariability.pdf', width=8, height=4 ,pointsize = 10)
par(mfrow=c(1,2));
par(mar=c(4, 4, 4, 4))
par(cex=1);
hist(unvariable, n=50, main="", xlab="Pearson's correlation between Col-0 and pif4-101")
plot(r_mat_pif42[1,], r_mat_pif42[2,], xlab="random network", ylab="predicted network", pch=20, col=rgb(0.5, 0.5, 0.5, 0.15), xlim=c(-0.8, 1), ylim=c(-0.8, 1))
#points(r_mat_pif4[1,], r_mat_pif4[2,], col="red", pch=19)
#points(r_mat_pif4[1,], r_mat_pif4[3,], col="cornflowerblue", pch=19)
points(r_mat_pif42[1,which(unvariable>0.97)], r_mat_pif42[2,which(unvariable>0.97)], col=rgb(1, 0, 0, 0.2), pch=20)
#points(r_mat_pif42[1,which(unvariable<0.2)], r_mat_pif42[2,which(unvariable<0.2)], col="blue")
abline(c(0,1))
abline(v=0) 
abline(h=0)
#legend(0.1, -0.3, c("randomly selected genes with G-boxes", "genes with PIF4 as predicted target", "PIF4 removed from predicted network"), pch=c(20, 19, 19), col=c("grey", "red", "cornflowerblue"), cex=0.8, bty='n')

dev.off()

