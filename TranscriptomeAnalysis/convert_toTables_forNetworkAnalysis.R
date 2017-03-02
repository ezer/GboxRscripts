#put batches in list
batches=list(batch1, batch2, batch3, batch4, batch5)

#gene names to extract-- binding sites
bindees=as.character(names_with_Gbox)
#bindees=as.character(names_with_HSE[,1])
#gene names to extract-- TFs
binders=as.character(read.table("../tf_names/gbox_TFs.txt")[,1])
#binders=as.character(read.table("../tf_names/hse_binders.txt")[,1])
allGenesToInclude=unique(c(binders, bindees))

#for each batch, extract gene names, merge values, make sure to remove missing binders
batches_zscore_binders_bindees_only=lapply(batches, function(batch){
  
  #isolate just relevant genes and put them in right order and assign proper rownames
  mats=lapply(batch, function(b){
    #isolate correct genes
    genes=as.character(b[,1])
    start=regexpr("AT", genes[1])
    gene_names_proper=substring(genes, start, start+9)
    b=b[which(gene_names_proper %in% allGenesToInclude),]
    
    #get genes in correct order
    genes=as.character(b[,1])
    start=regexpr("AT", genes[1])
    gene_names_proper=substring(genes, start, start+9)
    rownames(b)=gene_names_proper
    b=b[allGenesToInclude, ]
    b=b[which(as.character(rownames(b)) %in% allGenesToInclude),]
  })
  

   rowname_freqs=(table(unlist(lapply(mats, function(i){rownames(i)}))))
  
   allRowNamesToInclude=names(rowname_freqs)[which(rowname_freqs==length(mats))]
  
   mats=lapply(mats, function(b){
     b[allRowNamesToInclude,]
   })
   
   
   
   combined_mats=do.call("cbind", mats)
   rn=rownames(combined_mats)
   combined_mats=apply(combined_mats, 2, function(i){as.numeric(i)})
   rownames(combined_mats)=rn
   #get rid on non-numeric cols
   combined_mats=combined_mats[,which(apply(combined_mats,2, function(i){!anyNA(i)}))]
   
   #turn to z-score
  a= apply(combined_mats, 1, function(i){zscore(i)})
   rownames(a)=colnames(combined_mats)
   colnames(a)=rownames(combined_mats)
   a
  })

heatmap(batches_zscore_binders_bindees_only[[1]], Rowv=NA)

#save to file

for(i in c(1:length(batches_zscore_binders_bindees_only))){
  #table
  write.table(batches_zscore_binders_bindees_only[[i]], file=paste("batch_", i, "zscores.txt", sep=""), quote=FALSE)
  #table of all TF locations
  write.table(which(colnames(batches_zscore_binders_bindees_only[[i]]) %in% binders), file=paste("batch_", i, "tf_ids", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
  #names of all genes
   write.table(colnames(batches_zscore_binders_bindees_only[[i]]), quote=FALSE, row.names=FALSE, col.names = FALSE,file=paste("geneNames_batch_", i, "zscores.txt", sep=""))
}

# for(i in c(1:length(batches_zscore_binders_bindees_only))){
#   #table
#   write.table(batches_zscore_binders_bindees_only[[i]], file=paste("batch_", i, "zscores_hse.txt", sep=""), quote=FALSE)
#   #table of all TF locations
#   write.table(which(colnames(batches_zscore_binders_bindees_only[[i]]) %in% binders), file=paste("batch_", i, "tf_ids_hse", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
#   #names of all genes
#   write.table(colnames(batches_zscore_binders_bindees_only[[i]]), quote=FALSE, row.names=FALSE, col.names = FALSE,file=paste("geneNames_batch_", i, "zscores_hse.txt", sep=""))
# }



#check that it makes sense
plot(batches_zscore_binders_bindees_only[[2]][1:8,6])

id=which(as.character(batches[[2]][[1]][,1])==colnames(batches_zscore_binders_bindees_only[[2]])[6])

plot(as.numeric(batches[[2]][[1]][id,2:9]))

plot(as.numeric(batches[[2]][[1]][id,2:9]), batches_zscore_binders_bindees_only[[2]][1:8,6])

#names of TFs that are there