
# INPUT
# S : surrogate S
# Y : labels containing NA
# OUTPUT
# data.q : dataframe(S.new,Y) where S.new is the quantile transform of S
dat.quantile.transf = function(S,Y){
  f <- ecdf(S)
  S.new=f(S)
  dat.q = cbind(S.new,Y)
  dat.q
}

# INPUT
# S : surrogate S
# Y : labels containing NA
# OUTPUT
# data.new : dataframe(S.new,Y) where S.new  is obtained by taking bins of size h and 
# replacing the data S by its mean in the bin
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
