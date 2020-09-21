
#' Normalize the data using quantiles
#' @param S all surrogates S
#' @param Y labels containing NA
#' @return dataframe(S.new,Y) where S.new is the quantile transform of S
#' @export
dat.quantile.transf = function(S,Y){
  f <- ecdf(S)
  S.new=f(S)
  dat.q = cbind(S.new,Y)
  dat.q
}

#' Reduce the number of unique values of S. Take bins of size h and replacing the data S by its mean in the bin
#' @param S all surrogates S
#' @param Y labels containing NA
#' @return dataframe(S.new,Y) where S.new is the transform of S
#' @export
dat.new.FUN = function(S,Y,h){
  # h number of levels
  C = cut(S,h)
  dat.rank = data.frame(cbind(S,C))
  bin.means = aggregate(S ~ C,dat.rank,mean)
  length(unique(sort(dat.rank[,2])))
  # replace each rank value by respective mean
  r = as.factor(dat.rank[,2])
  m = bin.means[,2]
  S.new = m[r]
  dat.new = cbind.data.frame(S.new,Y)
  dat.new
}
