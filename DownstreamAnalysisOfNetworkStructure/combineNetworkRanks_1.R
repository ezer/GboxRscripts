
setwd("Data/networks")
pif4="AT2G43010"
CLR=read.table("gbox_minus_pif_CLR.txt")
gbox_genes=as.character(unique(CLR[,2]))
Tigress=read.table("gbox_minus_pif_tigress.txt", header=TRUE)
setwd("../..")
names_tigress=as.character(read.table("gbox_minus_pif_allrelatedgene_names.txt")[[1]])
setwd("Data/networks")
Tigress[,1]=names_tigress[as.numeric(Tigress[,1])]
Tigress[,2]=names_tigress[as.numeric(Tigress[,2])]
Genie=read.table("gbox_minus_pif_Genie3.txt")

names_genie=sapply(which(Genie[,3]!=0), function(i){
  paste(as.character(Genie[i,1]), as.character(Genie[i,2]))
})

names_clr=apply(CLR, 1, function(i){
  paste(as.character(i[1]), as.character(i[2]))
})

names_tigress=apply(Tigress, 1, function(i){
  paste(as.character(i[1]), as.character(i[2]))
})


order=sapply(names_clr, function(i){
  which(names_genie == i)
})

order_tig=sapply(names_tigress, function(i){
  which(names_genie == i)
})


plot(order, pch=19, col=rgb(0,0,1,0.005), xlab="CLR rank", ylab="Genie3 rank")
plot(order_tig, pch=19, col=rgb(0,0,1,0.005), xlab="Tigress rank", ylab="Genie3 rank")


unique_names=unique(c(names_genie, names_tigress, names_clr))

max_val=20000
avg_rank=sapply(unique_names, function(i){
  r1=which(names_genie[1:max_val]== i)
  if(length(r1)!=1){
    r1=max_val
  }
  
  r2=which(names_clr[1:max_val]== i)
  if(length(r2)!=1){
    r2=max_val
  }
  
  r3=which(names_tigress[1:max_val]== i)
  if(length(r3)!=1){
    r3=max_val
  }
  
  (r1+r2+r3)/3.0
})
hist(avg_rank[which(avg_rank<20000)], n=40, xlab="avg rank", ylab="frequency", main="")

save(avg_rank, file="avg_rank.RData")










