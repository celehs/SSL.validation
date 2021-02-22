

newfunction=function(roc,col.nm,value){
  cutoff.line = roc[which.min(abs(roc[,col.nm]-value)), ]
  cutoff = cutoff[cut]
}

#roc : num [1:99,1:6] roc table 
#col.nm : quantity of interest char
#value : real number
cutoff.choose=function(roc,col.nm,value){
  roc = data.frame(roc)
  cutoff.line = roc[which.min(abs(roc[,col.nm]-value))[1], ]
  cutoff = cutoff.line$cut
  res = list(cutoff,cutoff.line)
  res
}