setwd("Data/networks")
#save(allRNA, edges_networkSizes, file="testRobustness.RData")
load("avg_rank.RData")

######Step 1: draw clusters from Wild Type data (only show genes that are connected to the main network)

edges=sapply(c(1:length(avg_rank)), function(i){
  c(substr(names(avg_rank)[i], 1, 9), 
    substr(names(avg_rank)[i], 11, 19),
    avg_rank[i])
})
edges=t(edges)
edges[,3]=as.numeric(edges[,3])
######

#for each rna_seq_filtered, extract only the genes in network and put them in the same order
allGeneNames=unique(c(edges[,1], edges[,2]))
rna_seq_filtered_ordered=lapply(c(1:length(rna_seq_filtered)), function(i){
  print(i)
  if(substring(rna_seq_filtered[[i]][1,1], 1, 4)=="Gene"){
    b=substring(rna_seq_filtered[[i]][,1], 6, 14) 
  }else{
    b=rna_seq_filtered[[i]][,1]
  }
  a=rna_seq_filtered[[i]][which(as.character(b) %in% allGeneNames),]
  rownames(a)=b[which(as.character(b) %in% allGeneNames)]
  a=as.matrix(a[allGeneNames,2:length(a[1,])])#numeric bits only
  
  #remove columns that are just NAs
  
  a=a[,which(apply(a, 2, function(j){
    !(all(is.na(as.numeric(j))))
  }))]
  oldCols=colnames(a)
  a=apply(a, 1, function(j){
    zscore(as.numeric(j))}) #zscore for each row
  rownames(a)=oldCols
  a
})
names(rna_seq_filtered_ordered)=names(rna_seq_filtered)

allRNA=do.call("rbind", rna_seq_filtered_ordered)
#get 6 clusters
library(gplots)
a=heatmap.2(allRNA, trace="none")
hc <- as.hclust( a$colDendrogram )
a=cutree( hc, k=7 )
cols=c("cornflowerblue", "darkmagenta", "darkgreen", "grey", "darkblue", "lightgreen", "brown")
b=heatmap.2(allRNA, trace="none", ColSideColors = cols[a])
####draw this for WT genes only
allRNA_WT=do.call("rbind", lapply(c("col0_22", "col0_27", "ler0_22", "ler0_27"), function(i){rna_seq_filtered_ordered[[i]]}))

adjustColnames=substring(rownames(allRNA_WT), 5, nchar(rownames(allRNA_WT))-4)
col_temps=sapply(substring(rownames(allRNA_WT), nchar(rownames(allRNA_WT))-2, nchar(rownames(allRNA_WT))), function(i){if(i=="22c"){"blue"}else{"red"}})
sub1=c("ZT.4", "ZT.2", "ZT.12", "ZT.8")
sub2=c("ZT20",   "ZT22",   "ZT12", "ZT16")
for(i in c(1:4)){
  adjustColnames=gsub(sub1[i], sub2[i], adjustColnames)
}
rownames(allRNA_WT)=adjustColnames

cols=c("orange", "cornflowerblue", "tan", "darkgreen", "lightgreen", "red", "pink")
heatmap.2(allRNA, trace="none", ColSideColors = cols[a])
pdf(paste("Figure_clusterWT_allNetworkGenes.pdf", sep=""), width=7, height=7,pointsize = 10) 
par(mfrow=c(1,1));
par(mar=c(4, 4, 10, 10)+0.1);
par(cex=0.9)
c=heatmap.2(allRNA_WT[,b$colInd], trace="none", Colv=NA, dendrogram="none", Rowv=NA, margins=c(8,12), ColSideColors = cols[a[b$colInd]], labCol=NA, colRow=col_temps, add.expr=c(mtext("col-0", 4, 3.5, adj=0.75), mtext("ler", 4, 3.5, adj=0.25)))
mtext("col-0", 4, 1, adj=0.6)
mtext("ler", 4, 1, adj=0.1)

orderOfgenes=colnames(allRNA_WT)[b$colInd]
colorsOfGenes=cols[a[b$colInd]]


legend("bottomleft", legend=c("night", "late night+day", "night+early day", "early night", "day", "none", "dawn"), fill=cols[c(5, 7, 4, 2, 6, 1, 3)], bty="n", cex=0.7)

#legend(-200,1000, c("night", "late night+day", "night+early day", "early night", "day", "none", "dawn"), fill=cols[c(5, 7, 4, 2, 6, 1, 3)], bty="n", cex=0.7)
dev.off()





#######Step 2: draw full clusters, not just wild type

#########################################################################
#re-organise the stuff in combinedTable_norm to be more visually pleasing
############


orderOfSamples=c("col0_22",
                 "col0_27",
                 "ler0_22",
                 "ler0_27",
                 "elf3_22",
                 "elf3_27",
                 "lux_22",
                 "lux_27",
                 "arp6_22",
                 "arp6_27",
                 "dek3_17",
                 "dek3_27",
                 "DEKOX_17",
                 "DEKOX_27",
                 "pif4_22",
                 "pif4_27",
                 "phy_22",
                 "phy_27",
                 "hos1_22",
                 "hos1_27",
                 "ss4_22",
                 "ss4_27"
)

colCodes=c(137, 135, 375, 373, 150, 148, 145, 143, 125, 123, 132, 130, 30, 27, 258, 256, 614, 612, 99, 96, 588, 586)

#Get Z-score table
#extract genes from transcriptomes and generate giant single table
#transcriptimes=rna_seq_filtered
extractedSegment=sapply(orderOfSamples, function(m){
  mat=transcriptomes[[m]]
  #make sure gene names start with AT
  at=regexpr("AT", mat[1,1])
  genes=substring(mat[,1], at, at+9)
  
  #get it in the same order as the genes in box
  order=sapply(box, function(i){which(as.character(genes)==as.character(i))})
  print(m)
  print(dim(order))
  #extract that subsection
  if(m=="ss4_27" | m=="ss4_22"){
    a= mat[order,2:11] 
  }else{
    a=mat[order,2:length(mat[1,])]
  }
  rownames(a)=genes[order]
  a
})

combinedTable=do.call("cbind", extractedSegment)

#normalise table to Z-score
#rownames(combinedTable)=box
#combinedTable=combinedTable[,2:length(combinedTable[1,])]
combinedTable_numeric=combinedTable[,which(!apply(combinedTable, 2, function(i){anyNA(as.numeric(i))}))]
combinedTable_norm=apply(as.matrix(combinedTable_numeric), 1, function(i){ zscore(as.numeric(i))})

cols=colors()[colCodes]
repCount=sapply(extractedSegment, function(i){length(which(sapply(i[1,], function(j){is.numeric(j)})))})
colsCol=unlist(sapply(c(1:length(extractedSegment)), function(i){rep(cols[i], repCount[i])}))
key=unlist(sapply(c(1:length(extractedSegment)), function(i){rep(names(repCount)[i], repCount[i])}))
pdf('Figure_heatmapsCluster_rearranged_byType_v3.pdf', width=10, height=10,pointsize = 10) 
par(mfrow=c(1,1));


heatmap(combinedTable_norm[,orderOfgenes], Rowv=NA, Colv=NA, RowSideColors = colsCol, ColSideColors = colorsOfGenes,  labRow = NA, labCol = NA)

#heatmap(combinedTable_norm[which(key!="ss4_17_" & key!="col17_"),], Rowv=NA, RowSideColors = colsCol[which(key!="ss4_17_" & key!="col17_")])


dev.off()