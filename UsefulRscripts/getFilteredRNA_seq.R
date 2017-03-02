getFilteredRNA_seq <-function(batchName, experiments, filtered=TRUE){


setwd(batchName)

#read in each file
file_names=paste(experiments, "_tpm", sep="")
raw_data=lapply(file_names, function(i){
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


if(filtered){
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
setwd("..")
rna_seq_filtered
}else{
  setwd("..")
  rna_seq
}
}
