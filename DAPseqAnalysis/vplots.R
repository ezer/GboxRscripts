setwd("DAPseqAnalysis")

library("Rsamtools")
library("GenomicAlignments")
#load MNase data
files=c("Input_HTA11_17C_T0-38167253_raw_trimmo_paired_truseq3-PE-2_2_10_5_1_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam",
        "Input_HTA11_17C_T8-38193198_raw_trimmo_paired_truseq3-PE-2_2_10_5_1_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam")

#gene_coords <- GRanges (seqnames = chr, ranges = IRanges(start = cds_start, end = cds_end))

bam_fields <- c("rname", "pos", "qwidth", "mpos", "isize") #I believe these are the minimal parameters required

param <- ScanBamParam(what = bam_fields)

#bam_file <- "2_33C_03_chip-seq_20141101_col_27c15min_dark_mnase_raw_trimmo_paired_truseq3-PE_2_10_5_1_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam"
bams=lapply(files, function(bam_file){
  readGAlignments(bam_file, param = param)
})

vplot_components=lapply(c(1:2), function(f){
  lapply(c(1:5), function(chr){
    arr=bams[[f]]@elementMetadata
    chrRows=which(as.character(arr[,"rname"])==as.character(chr))
    pos_centre=arr[chrRows,"pos"]+0.5*arr[chrRows, "isize"]
    fragLen=arr[chrRows, "isize"]
    cbind(pos_centre, fragLen)
  })
})

######### make *binding* clusters, rather than expression clusters
cls=lapply(c(1:4), function(i){
  as.character(read.table(paste("cl", i, ".txt", sep=""))[,1])
})

all_locations_by_binding=lapply(c(1:4), function(cluster){
  sapply(cls[[cluster]], function(i){
    #get info
    metadata=gbox_geneTable[which(as.character(gbox_geneTable[,1])==as.character(i)),2]
    chr=substring(metadata, 3, 6)
    coord=as.numeric(strsplit(substring(metadata, 8, nchar(as.character(metadata))-19), "-")[[1]])
    
    #get coordinates of perfect G-boxes within sequence
    gboxSplit=as.numeric(gregexpr("CACGTG", as.character(gbox_geneTable[which(as.character(gbox_geneTable[,1])==as.character(i)),3]))[[1]])
    
    #get forward or reverse
    if(length(grep("FORWARD", as.character(metadata)))>0){
      c(chr, coord[1]+gboxSplit-1)  
    }else{
      c(chr, coord[2]-gboxSplit+1)  
    }
    
    #
  })
})


pdf(paste("Figure_Gbox_vplots_byBindingPattern.pdf", sep=""), width=10, height=3.5,pointsize = 10) 
par(mfrow=c(1,3));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

for(cl in c(1:2)){
  plot(c(-flank, flank), c(50, 400), col="white", main=paste("binding pattern", cl), xlab="centre of fragment relatve to G-box", ylab="fragment size")
  for(d in all_locations_by_binding[[cl]]){ 
    chr=as.numeric(substr(d[1], 4, 4))
    for(coord in as.numeric(d[2:length(d)])){
      ids=which(vplot_components[[f]][[chr]][, 1]>(coord-flank) & vplot_components[[f]][[chr]][, 1]<(coord+flank))
      points(vplot_components[[f]][[chr]][ids, 1]-coord,
             vplot_components[[f]][[chr]][ids, 2], pch=20, col=rgb(0, 0, 1, 0.01))
          }
  }
 }

dev.off()

