
setwd("../genes")

#read in all of the gene sequences
geneSeqs=lapply(c(1:10), function(i){
  read.table(paste("genes_", i, sep=""), quote="", comment.char="", header=FALSE, sep="\t")
})

geneSeqs_unlisted=do.call("rbind", geneSeqs)
#save dataset
save(geneSeqs_unlisted, file="geneSeqs_unlisted.RData")

#####################################################
######################Find genes with perfect G-boxes
###################### CACGTG
#find all sequences with perfect G-boxes
locationsOfGbox=regexpr("CACGTG", as.character(geneSeqs_unlisted[,3]))
#Draw a histogram of the distribution of the G-box locations
allGboxes=unlist(gregexpr("CACGTG", as.character(geneSeqs_unlisted[,3])))
hist(1000-allGboxes[which(allGboxes>0)], xlab="bp from TSS", ylab="frequency of G-box", main="", n=60)


#This analysis suggests that there is an enrichment of G-boxes within 500bp of the TSS
#get all the 500bp sequences
geneSeqs_unlisted_500bp=geneSeqs_unlisted
geneSeqs_unlisted_500bp[,3]=substring(as.character(geneSeqs_unlisted_500bp[,3]), 501, 1000)

#select genes with perfect G-boxes for future analysis
allGboxes=unlist(gregexpr("CACGTG", as.character(geneSeqs_unlisted_500bp[,3])))
hist(500-allGboxes[which(allGboxes>0)], xlab="bp from TSS", ylab="frequency of G-box", main="",n=30)

#save gene list of genes that have perfect G-boxes
names_with_Gbox=geneSeqs_unlisted_500bp[
  which(regexpr("CACGTG", as.character(geneSeqs_unlisted_500bp[,3]))>0),1]

write.table(names_with_Gbox, quote=FALSE, row.names = FALSE, col.names = FALSE, file="genesWithGboxes.txt")



#####################################################
######################Find genes with perfect W-boxes
###################### TTGAC(C/T)

#Draw a histogram of the distribution of the W-box locations
allWboxes=c(unlist(gregexpr("TTGACC", as.character(geneSeqs_unlisted[,3]))), 
            unlist(gregexpr("TTGACT", as.character(geneSeqs_unlisted[,3]))),
            unlist(gregexpr("AGACAA", as.character(geneSeqs_unlisted[,3]))),
            unlist(gregexpr("GGACAA", as.character(geneSeqs_unlisted[,3]))))
hist(1000-allWboxes[which(allWboxes>0)], xlab="bp from TSS", ylab="frequency of G-box", main="")

allWboxes=c(unlist(gregexpr("TTGACC", as.character(geneSeqs_unlisted_500bp[,3]))), 
            unlist(gregexpr("TTGACT", as.character(geneSeqs_unlisted_500bp[,3]))),
            unlist(gregexpr("AGACAA", as.character(geneSeqs_unlisted_500bp[,3]))),
            unlist(gregexpr("GGACAA", as.character(geneSeqs_unlisted_500bp[,3]))))
hist(500-allWboxes[which(allWboxes>0)], xlab="bp from TSS", ylab="frequency of G-box", main="")


names_with_Wbox_ids=which(sapply(as.character(geneSeqs_unlisted_500bp[,3]), function(i){
  a=sapply(c("TTGACC", "TTGACT", "AGACAA", "GGACAA"), function(j){
    regexpr(j, i)>0
  })
  length(which(a))>0
}))

names_with_Wbox=as.character(geneSeqs_unlisted_500bp[names_with_Wbox_ids,1])

write.table(names_with_Wbox, quote=FALSE, row.names = FALSE, col.names = FALSE, file="genesWithWboxes.txt")

#####################################################
######################Find genes with perfect HSE
###################### this time, I'm just going to use all 1000bp, because there is no <500bp bias in location of sequence
######################  I am also allowing for 1 mistake, because there are only <300 examples of perfect TTC..GAA..TTC and GAA..TTC..GAA


#maybe allow for one mistake?
motif_core=c("TTC..GAA..TTC", "GAA..TTC..GAA")
motifs=unlist(lapply(motif_core, function(i){
  sapply(c(1:3, 6:8, 11:13), function(j){
    paste(substring(i, 1, j-1), ".", substring(i, j+1,nchar(i)), sep="")
  })
}))

names_with_HSE_ids=which(sapply(as.character(geneSeqs_unlisted[,3]), function(i){
  #sapply(c("TTC..GAA..TTC", "GAA..TTC..GAA"), function(j){
  #a=sapply(c("GAA..TTC", "TTC..GAA"), function(j){
  a=sapply(motifs, function(j){ 
    regexpr(j, i)>0
  })
  length(which(a))>0
}))

names_with_HSE=as.character(geneSeqs_unlisted[names_with_HSE_ids,1])

write.table(names_with_HSE, quote=FALSE, row.names = FALSE, col.names = FALSE, file="genesWithHSE.txt")

names_with_HSE=read.table("genesWithHSE.txt")






