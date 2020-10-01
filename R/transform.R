
dat.quantile.transf = function(S,Y){
  f <- ecdf(S)
  S.new=f(S)
  dat.q = cbind(S.new,Y)
  dat.q
}

#' Choose the cutoff corresponding to desired 
#' @param roc roc table, list 
#' @param col.nm column name of the roc table, string
#' @param value value of desired 
#' @return list with cutoff value and line of the roc table for this cutoff
#' @export
cutoff.choose=function(roc,col.nm,value){
  roc = data.frame(roc)
  cutoff.line = roc[which.min(abs(roc[,col.nm]-value))[1], ]
  cutoff = cutoff.line$cut
  res = list(cutoff,cutoff.line)
  res
}

dat.new.FUN = function(S,Y,h){
  C = cut(S,h)
  dat.rank = data.frame(cbind(S,C))
  bin.means = aggregate(S ~ C,dat.rank,mean)
  length(unique(sort(dat.rank[,2])))
  r = as.factor(dat.rank[,2])
  m = bin.means[,2]
  S.new = m[r]
  dat.new = cbind.data.frame(S.new,Y)
  dat.new
}
