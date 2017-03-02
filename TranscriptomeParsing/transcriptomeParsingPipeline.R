#######################
# Before identifying subset of genes, convert fpkms into tpms
#######################

#setwd("~/Documents/GboxWebsite2/TranscriptomeParsing")
setwd("~/Documents/GboxWebsite2/Data")

#turn lux into tpms
#lux=read.table("textfiles/lux_fpkm", quote="", header=TRUE, comment.char="", sep="\t")
#lux_asTPM=cbind(lux[,1:3], fpkm_to_tpm(lux[, 4:length(colnames(lux))]))
#colnames(lux_asTPM)=sub("FPKM", "TPM", colnames(lux_asTPM))
#write.table(lux_asTPM, file="textfiles/lux_tpm", sep="\t", quote=FALSE)

#turn phya-e into tpms
#phy=read.table("textfiles/phya-e_fpkm", quote="", header=TRUE, comment.char="", sep="\t")
#phy_asTPM=cbind(phy[,1:3], fpkm_to_tpm(phy[, 4:length(colnames(phy))]))
#colnames(phy_asTPM)=sub("FPKM", "TPM", colnames(phy_asTPM))
#write.table(phy_asTPM, file="textfiles/phy_tpm", sep="\t", quote=FALSE)


#######################
# Use only data batch 1 to determine the genes that should be under study
#######################
setwd("textfiles")
setwd("batch1_20141105_20140901")

#read in each file
experiments=c("col0", "ler0", "elf3", "lux", "phy")
file_names=paste(experiments, "_tpm", sep="")
raw_data=sapply(file_names, function(i){
  read.table(i, quote="", header=TRUE, comment.char="", sep="\t")
})
names(raw_data)=experiments

#separate out by temperature
rna_seq=list()
for(mat_name in names(raw_data)){
  unique_temps=unique(getTemperature_2digits(colnames(raw_data[[mat_name]])[4:length(colnames(raw_data[[mat_name]]))]))
  temps=getTemperature_2digits(colnames(raw_data[[mat_name]])[1:length(colnames(raw_data[[mat_name]]))])
  for(t in unique_temps){
    rna_seq[[paste(mat_name, "_", t, sep="")]]=cbind(raw_data[[mat_name]][,1], raw_data[[mat_name]][,which(temps==t)])
  }
}

#remove columns that have no data (NA or Div/0)
#remove rows where gene is <1 TPM more than half the time
rna_seq_filtered=lapply(rna_seq, function(mat){
  temp=cbind(mat[,1], mat[,1+which(apply(as.matrix(mat[,2:length(colnames(mat))]), 2, function(i){
    a=as.numeric(i)
    length(which(!is.na(a) & !is.nan(a) & !is.infinite(a)))>1
  }))])
  
  temp[which(apply(temp, 1, function(i){
    a=as.numeric(i)
    length(which(a>1))>(length(a)/2)
  })),
  ]
  
  
})
names(rna_seq_filtered)=names(rna_seq)

#find set of genes that are found in every data set
#filter for genes that have >1 TPM in at least half the samples in each dataset 
allGeneNames=unlist(sapply(rna_seq_filtered, function(i){
  names=as.character(i[,1])
  ind=regexpr("AT", names)[1]
  unique(substr(names, ind, ind+9))
  
}))
freqNames=table(allGeneNames)
highEnoughExpressionNames=names(freqNames)[which(freqNames==length(rna_seq_filtered))]

#save lists of these genes for sequence analysis
#save these names in groups of 1,900 because the maximum number of downloads for TAIR10 is 2,000
ranges=c(0, c(1:8)*1900, length(highEnoughExpressionNames))
allGenesListed=c()
for(i in c(2:length(ranges))){
  write.table(highEnoughExpressionNames[c((1+ranges[i-1]):ranges[i])], quote=FALSE, row.names = FALSE, col.names = FALSE, file=paste("geneNames_", (i-1), sep=""))
  allGenesListed=c(allGenesListed, highEnoughExpressionNames[c((1+ranges[i-1]):ranges[i])])
}


#create files for each of the other "batches"
setwd("..")
batch1=getFilteredRNA_seq("batch1_20141105_20140901", c("col0", "ler0", "elf3", "lux", "phy"))

batch2=getFilteredRNA_seq("batch2_20141101", c("hos1", "arp6", "pif4"))

batch3=getFilteredRNA_seq("batch3_20160114", c("col0ForDek", "dek3", "DEKOX"))

batch4=getFilteredRNA_seq("batch4_17C", c("col0_forss4_17", "ss4_17"))

batch5=getFilteredRNA_seq("batch5_ss4", c("col0_forss4", "ss4"))

batch_labels=list()
batch_labels[[1]]=c("col0", "ler0", "elf3", "lux", "phy")
batch_labels[[2]]=c("hos1", "arp6", "pif4")
batch_labels[[3]]=c("col0ForDek", "dek3", "DEKOX")
batch_labels[[4]]=c("col0_forss4_17", "ss4_17")
batch_labels[[5]]=c("col0_forss4", "ss4")
rna_seq_filtered_2=c(batch1, batch2, batch3, batch4, batch5)






















#########old version
#read in each file
experiments=c("arp6", "hos1", "pif4", "dek3", "DEKOX", "col0", "elf3", "ler0", "lux", "phya-e", "col17", "ss4_17", "ss4", "yhb27", "yhb22")
file_names=paste(experiments, "_tpm", sep="")
raw_data=sapply(file_names, function(i){
  read.table(paste("textfiles/", i, sep=""), quote="", header=TRUE, comment.char="", sep="\t")
})
names(raw_data)=experiments

#separate out by temperature
rna_seq=list()
for(mat_name in names(raw_data)){
  unique_temps=unique(getTemperature_2digits(colnames(raw_data[[mat_name]])[4:length(colnames(raw_data[[mat_name]]))]))
  temps=getTemperature_2digits(colnames(raw_data[[mat_name]])[1:length(colnames(raw_data[[mat_name]]))])
  for(t in unique_temps){
    rna_seq[[paste(mat_name, "_", t, sep="")]]=cbind(raw_data[[mat_name]][,1], raw_data[[mat_name]][,which(temps==t)])
  }
}

  
#remove columns that have no data (NA or Div/0)
#remove rows where gene is <1 TPM more than half the time
rna_seq_filtered=lapply(rna_seq, function(mat){
  temp=cbind(mat[,1], mat[,1+which(apply(as.matrix(mat[,2:length(colnames(mat))]), 2, function(i){
    a=as.numeric(i)
    length(which(!is.na(a) & !is.nan(a) & !is.infinite(a)))>1
}))])
  
  temp[which(apply(temp, 1, function(i){
    a=as.numeric(i)
    length(which(a>1))>(length(a)/2)
  })),
       ]
  
  
  })
names(rna_seq_filtered)=names(rna_seq)

#find set of genes that are found in every data set
#filter for genes that have >1 TPM in at least half the samples in each dataset 
allGeneNames=unlist(sapply(rna_seq_filtered, function(i){
  names=as.character(i[,1])
  ind=regexpr("AT", names)[1]
  unique(substr(names, ind, ind+9))
  
}))
freqNames=table(allGeneNames)
highEnoughExpressionNames=names(freqNames)[which(freqNames==length(rna_seq_filtered))]


#save lists of these genes for sequence analysis
#save these names in groups of 1,900 because the maximum number of downloads for TAIR10 is 2,000
ranges=c(0, c(1:8)*1900, length(highEnoughExpressionNames))
allGenesListed=c()
for(i in c(2:length(ranges))){
  write.table(highEnoughExpressionNames[c((1+ranges[i-1]):ranges[i])], quote=FALSE, row.names = FALSE, col.names = FALSE, file=paste("geneNames_", (i-1), sep=""))
  allGenesListed=c(allGenesListed, highEnoughExpressionNames[c((1+ranges[i-1]):ranges[i])])
  }


#save transcriptome data
save(rna_seq_filtered, file="rna_seq_filtered.RData")



#######################
# After identifying subset of genes
######################

#filter out genes

#make giant table

#generate Z-scores

#make alternative tables

#cluster genes

