setwd("DAPseqAnalysis")


#######Look for overlaps of all perfect G-box regions
gbox_tf_overlaps=read.table("DAPseq_match_freq_gbox.txt") 
possibilities=gbox_tf_overlaps[which(gbox_tf_overlaps[,2]>300),]
possibilities[order(possibilities[,2], decreasing=TRUE),]

gbox_regions=read.table("perfectGboxRegions.bed", strip.white=TRUE)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
gbox_names=apply(gbox_regions, 1, function(i){paste(trim(i[1]), trim(i[2]), trim(i[3]), sep="_")})
tf_names=possibilities[,1]
freq_mat=matrix(0, nrow=length(tf_names), ncol=length(gbox_names))
rownames(freq_mat)=tf_names
colnames(freq_mat)=gbox_names

#load up other data file...
con <- file("DAPseq_match_gbox.txt") 
open(con);
#results.list <- list();
current.line <- 1
header=""
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  #print(as.character(line))
  if(substring(as.character(line), 1, 4)=="proc"){
    #print("proc")
    header=substring(as.character(line), 60, nchar(as.character(line))-28)
    print(header)
  }else{
    if(header %in% as.character(tf_names)){
      print(line)
      i=strsplit(as.character(line), "\t")
      print(i)
      gbox_id=paste(i[[1]][1], i[[1]][2], i[[1]][3], sep="_")
      print(gbox_id)
      freq_mat[header, gbox_id]=1
    }
  }
  
  
  
  # results.list[[current.line]] <- as.integer(unlist(strsplit(line, split=" ")))
  current.line <- current.line + 1
} 
close(con)



cols=c("cornflowerblue", "darkmagenta", "red", "orange", "yellow", "green", "darkblue", "darkgreen", "pink", "brown")
names(cols)=c("bHLH", "bZIP", "C2C2dof", "HB", "ZFHD", "MYB", "BZR", "NAC", "AP2EREBP", "WRKY")
col_col=sapply(rownames(freq_mat), function(i){
  i_spl=strsplit(i, "/")[[1]][2]
  foundMatch=which(sapply(names(cols), function(j){length(grep(j, i_spl))==1}))
  if(length(foundMatch)!=1){
    "grey"
  }else{
    cols[foundMatch]
  }
})
pdf('Figure_clusterGboxBinders.pdf', width=8, height=8 ,pointsize = 10)
par(mfrow=c(1,1));
par(mar=c(4, 4, 4, 4))
par(cex=0.2);
heatmap(-freq_mat, scale="none", labRow=NA, labCol=NA, xlab = "G-box containing promoter region", ylab="TF overlapping with >15% of regions", RowSideColors = col_col)
dev.off()

########################The same as above but have ALL TFs


#######Look for overlaps of all perfect G-box regions
gbox_tf_overlaps=read.table("DAPseq_match_freq_gbox.txt") 
possibilities=gbox_tf_overlaps[which(gbox_tf_overlaps[,2]>0),]
possibilities[order(possibilities[,2], decreasing=TRUE),]

gbox_regions=read.table("perfectGboxRegions.bed", strip.white=TRUE)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
gbox_names=apply(gbox_regions, 1, function(i){paste(trim(i[1]), trim(i[2]), trim(i[3]), sep="_")})
tf_names=possibilities[,1]
freq_mat=matrix(0, nrow=length(tf_names), ncol=length(gbox_names))
rownames(freq_mat)=tf_names
colnames(freq_mat)=gbox_names

#load up other data file...
con <- file("DAPseq_match_gbox.txt") 
open(con);
#results.list <- list();
current.line <- 1
header=""
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  #print(as.character(line))
  if(substring(as.character(line), 1, 4)=="proc"){
    #print("proc")
    header=substring(as.character(line), 60, nchar(as.character(line))-28)
    print(header)
  }else{
    if(header %in% as.character(tf_names)){
      print(line)
      i=strsplit(as.character(line), "\t")
      print(i)
      gbox_id=paste(i[[1]][1], i[[1]][2], i[[1]][3], sep="_")
      print(gbox_id)
      freq_mat[header, gbox_id]=1
    }
  }
  
  
  
  # results.list[[current.line]] <- as.integer(unlist(strsplit(line, split=" ")))
  current.line <- current.line + 1
} 
close(con)



cols=c("cornflowerblue", "darkmagenta", "red", "orange", "yellow", "green", "darkblue", "darkgreen", "pink", "brown")
names(cols)=c("bHLH", "bZIP", "C2C2dof", "HB", "ZFHD", "MYB", "BZR", "NAC", "AP2EREBP", "WRKY")
col_col=sapply(rownames(freq_mat), function(i){
  i_spl=strsplit(i, "/")[[1]][2]
  foundMatch=which(sapply(names(cols), function(j){length(grep(j, i_spl))==1}))
  if(length(foundMatch)!=1){
    "grey"
  }else{
    cols[foundMatch]
  }
})


write.table(freq_mat, file="DAPseqResults_Gboxes.txt", quote=F)





pdf('Figure_clusterGboxBinders_legend.pdf', width=8, height=8 ,pointsize = 10)
par(mfrow=c(1,1));
par(mar=c(4, 4, 4, 4))
par(cex=1);
plot(c(0,1), c(0,1), col="white")
legend(0.5, 0.8, c(names(cols), "other"), fill=c(cols, "grey"), bty="n")
dev.off()

#save(freq_mat, cols, file="clusterGboxBinders.RData")


#now I can redo this, but just with all the bZIPs and bHLHs

gbox_tf_overlaps=read.table("DAPseq_match_freq_gbox.txt") 
gboxers=c("bHLH", "bZIP", "BZR", "lowQuality")
possibilities=gbox_tf_overlaps[which(sapply(as.character(gbox_tf_overlaps[,1]), function(i){
  i_spl=strsplit(i, "/")[[1]][2]
  foundMatch=which(sapply(gboxers, function(j){length(grep(j, i_spl))==1}))
  length(foundMatch)==1})),]

possibilities[order(possibilities[,2], decreasing=TRUE),]

gbox_regions=read.table("perfectGboxRegions.bed", strip.white=TRUE)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
gbox_names=apply(gbox_regions, 1, function(i){paste(trim(i[1]), trim(i[2]), trim(i[3]), sep="_")})
tf_names=possibilities[,1]
freq_mat=matrix(0, nrow=length(tf_names), ncol=length(gbox_names))
rownames(freq_mat)=tf_names
colnames(freq_mat)=gbox_names

#load up other data file...
con <- file("DAPseq_match_gbox.txt") 
open(con);
#results.list <- list();
current.line <- 1
header=""
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  #print(as.character(line))
  if(substring(as.character(line), 1, 4)=="proc"){
    #print("proc")
    header=substring(as.character(line), 60, nchar(as.character(line))-28)
    print(header)
  }else{
    if(header %in% as.character(tf_names)){
      print(line)
      i=strsplit(as.character(line), "\t")
      print(i)
      gbox_id=paste(i[[1]][1], i[[1]][2], i[[1]][3], sep="_")
      print(gbox_id)
      freq_mat[header, gbox_id]=1
    }
  }
  
  
  
  # results.list[[current.line]] <- as.integer(unlist(strsplit(line, split=" ")))
  current.line <- current.line + 1
} 
close(con)

cols=c("cornflowerblue", "darkmagenta", "darkblue", "orange")
names(cols)=c("bHLH", "bZIP", "BZR", "lowQuality")
col_col=sapply(rownames(freq_mat), function(i){
  i_spl=strsplit(i, "/")[[1]][2]
  foundMatch=which(sapply(names(cols), function(j){length(grep(j, i_spl))==1}))
  if(length(foundMatch)!=1){
    "grey"
  }else{
    cols[foundMatch]
  }
})
a=heatmap(-freq_mat, scale="none", labRow=NA, labCol=NA, xlab = "G-box containing promoter region", ylab="TF from G-box binding family", RowSideColors = col_col)
a=cutree(hclust(dist(t(freq_mat))), k=4)
pdf('Figure_clusterGboxBinders_bHLH_bZIP_BZR2_v2.pdf', width=8, height=8 ,pointsize = 10)
par(mfrow=c(1,1));
par(mar=c(4, 4, 4, 4))
par(cex=0.2);
v=heatmap(-freq_mat, scale="none", labRow=NA, labCol=NA, xlab = "G-box containing promoter region", ylab="TF from G-box binding family", ColSideColors=cols[a], RowSideColors = col_col)
dev.off()

pdf('Figure_clusterGboxBinders__bHLH_bZIP_BZR_legend.pdf', width=8, height=8 ,pointsize = 10)
par(mfrow=c(1,1));
par(mar=c(4, 4, 4, 4))
par(cex=1);
plot(c(0,1), c(0,1), col="white")
legend(0.5, 0.8, c(names(cols)[1:3], "low peak count bHLH"), fill=c(cols), bty="n")
dev.off()

save(freq_mat, cols, file="clusterGboxBinders_bHLH_bZIP_BZR.RData")







