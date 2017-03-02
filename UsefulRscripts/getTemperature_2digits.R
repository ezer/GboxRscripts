getTemperature_2digits <- function(names){
 
  a=regexpr("[123456789][1234567890][C]", toupper(names))
  substr(names,a, a+1)
  
}