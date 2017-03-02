fpkm_to_tpm <- function(mat){
  apply(mat, 2, function(i){
    i/(sum(i)/(10^6));
  })
}