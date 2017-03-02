zscore<-function(r){
  return ((r-mean(r, na.rm=TRUE))/sd(r, na.rm=TRUE))
}