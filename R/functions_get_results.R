analysis_run = function(Sv, St, Yt, Wt = NULL, b = 200){
  ## St: score in the labeled set
  ## Yt: outcome
  ## Sv: score in the unlabeled set
  ## Wt: optional vector of weights
  Y=c(Yt,rep(NA,length(Sv)))
  S=c(St,Sv)
  n.t = length(Yt); if(is.null(Wt)){Wt = rep(1, n.t)}
  ## get the estimates for entire data
  roc.sl = roc.superv(St,Yt)
  roc.sl$roc = unlist(cutoff.choose(roc.sl$roc,"fpr",0.05)[2])
  roc.ssl = roc.semi.superv(S,Y)
  roc.ssl$roc = unlist(cutoff.choose(roc.ssl$roc,"fpr",0.05)[2])
  set.seed(1234)
  ## perturbation 
  params = get_beta_params(0.25, 1/16)
  inds.resamp = sapply(1:b, function(jj) 4*rbeta(n.t, params$alpha, params$beta))
  
  n = nrow(inds.resamp)
  roc.ssl.pert = lapply(1:b, function(jj) tryCatch(roc.semi.superv.pert(S,Y, Wt*inds.resamp[,jj]), 
                                                   error=function(e) NA))
  for (i in 1:b){roc.ssl.pert[[i]]$roc = tryCatch(unlist(cutoff.choose(roc.ssl.pert[[i]]$roc,"fpr",0.05)[2]), error=function(e) NA)}
  
  roc.sl.pert = lapply(1:b, function(jj) tryCatch(roc.superv(St, Yt, wgti=Wt*inds.resamp[,jj]), 
                                                  error=function(e) NA))
  for (i in 1:b){roc.sl.pert[[i]]$roc = tryCatch(unlist(cutoff.choose(roc.sl.pert[[i]]$roc,"fpr",0.05)[2]), error=function(e) NA)}
  
  return(list(roc.sl = roc.sl, roc.ssl = roc.ssl, roc.ssl.pert = roc.ssl.pert, roc.sl.pert =roc.sl.pert))
}



analysis_summary = function(roc.sl, roc.ssl, roc.ssl.pert, roc.sl.pert){
  level=0.95
  b=length(roc.ssl.pert)
  
  auc.sl = roc.sl$auc
  roc.table.sl= roc.sl$roc
  auc.ssl = roc.ssl$auc
  roc.table.ssl = roc.ssl$roc
  
  auc.sl.list=NULL
  auc.ssl.list = NULL
  roc.sl.list = NULL
  roc.ssl.list = NULL
  for (i in 1:b){
    auc.sl.list = rbind(auc.sl.list, roc.sl.pert[[i]]$auc)
    auc.ssl.list = rbind(auc.ssl.list,roc.ssl.pert[[i]]$auc)
    roc.sl.list = rbind(roc.sl.list,roc.sl.pert[[i]]$roc)
    roc.ssl.list = rbind(roc.ssl.list,roc.ssl.pert[[i]]$roc)
  }
  
  sl.tab = cbind( rbind(auc.sl, sd(auc.sl.list, na.rm=T)), 
                  rbind(roc.table.sl, apply(roc.sl.list[complete.cases(roc.ssl.list),], 2, sd)))
  rownames(sl.tab)=c('Pt.est','sd')
  ssl.tab = cbind( rbind(auc.ssl, sd(auc.ssl.list, na.rm=T)), 
                   rbind(roc.table.ssl, apply(roc.ssl.list[complete.cases(roc.ssl.list),], 2, sd)))
  rownames(ssl.tab)=c('Pt.est','sd')
  
  return(list(sl.tab=sl.tab, ssl.tab=ssl.tab))
}




analysis_run_boot = function(Sv, St, Yt, Yv.na, Yv, Wt = NULL, b = 200){
  ## St: score in the labeled set
  ## Yt: outcome
  ## Sv: score in the unlabeled set
  ## Wt: optional vector of weights
  ## b: number of resamples
  ## ex. for subsampling set pert = FALSE and prop = 0.6 to subsample 60% of data wo replacement
  Y=Yv.na
  S=Sv
  n.t = length(Yt); if(is.null(Wt)){Wt = rep(1, n.t)}
  ## get the estimates for entire data
  roc.sl = roc.superv(St, Yt)
  roc.sl$roc = unlist(cutoff.choose(roc.sl$roc,"fpr",0.05)[2])
  
  roc.ssl = roc.semi.superv(S, Y)
  roc.ssl$roc = unlist(cutoff.choose(roc.ssl$roc,"fpr",0.05)[2])
  
  roc.c = roc.SSL.combine(Sv, Yv.na)
  roc.c$roc = unlist(cutoff.choose(roc.c$roc,"fpr",0.05)[2])
  
  roc.true = roc.superv(Sv, Yv)
  roc.true$roc = unlist(cutoff.choose(roc.true$roc,"fpr",0.05)[2])
  
  id.v = which(!is.na(Yv.na))
  
  ## bootstrap 
  yes.replace = TRUE
  inds.resamp = sapply(1:b, function(jj) sample(1:n.t, size = floor(1*n.t), replace = yes.replace))
  
  n = nrow(inds.resamp);
  roc.ssl.pert = lapply(1:b, function(jj) {inds.tmp = inds.resamp[,jj];
  tryCatch(roc.semi.superv(c(S[id.v[inds.tmp]], S[-id.v]), c(Y[id.v[inds.tmp]], Y[-id.v]) ), error=function(e) NA)})
  for (i in 1:b){roc.ssl.pert[[i]]$roc = tryCatch(unlist(cutoff.choose(roc.ssl.pert[[i]]$roc,"fpr",0.05)[2]), error=function(e) NA)}
  
  roc.c.pert = lapply(1:b, function(jj) {inds.tmp = inds.resamp[,jj];
  tryCatch(roc.SSL.combine(c(S[id.v[inds.tmp]], S[-id.v]), c(Y[id.v[inds.tmp]], Y[-id.v]) ), error=function(e) NA)})
  for (i in 1:b){roc.c.pert[[i]]$roc = tryCatch(unlist(cutoff.choose(roc.c.pert[[i]]$roc,"fpr",0.05)[2]), error=function(e) NA)}
  
  roc.sl.pert = lapply(1:b, function(jj) {inds.tmp = inds.resamp[,jj];
  tryCatch(roc.superv(St[inds.tmp], Yt[inds.tmp]), error=function(e) NA)})
  for (i in 1:b){roc.sl.pert[[i]]$roc = tryCatch(unlist(cutoff.choose(roc.sl.pert[[i]]$roc,"fpr",0.05)[2]), error=function(e) NA)}
  
  return(list(roc.sl = roc.sl, roc.ssl = roc.ssl, roc.c=roc.c, roc.true = roc.true, roc.sl.pert = roc.sl.pert, roc.ssl.pert =roc.ssl.pert, roc.c.pert=roc.c.pert))
}



analysis_summary_boot = function(roc.sl, roc.ssl, roc.c, roc.true, roc.sl.pert, roc.ssl.pert, roc.c.pert){
  b = length(roc.ssl.pert)

  auc.sl = roc.sl$auc
  roc.table.sl= roc.sl$roc
  auc.ssl = roc.ssl$auc
  roc.table.ssl = roc.ssl$roc
  auc.c = roc.c$auc
  roc.table.c = roc.c$roc
  auc.true = roc.true$auc
  roc.table.true = roc.true$roc
  
  auc.sl.list=NULL
  auc.ssl.list = NULL
  auc.c.list = NULL
  roc.sl.list = NULL
  roc.ssl.list = NULL
  roc.c.list = NULL
  for (i in 1:b){
    auc.sl.list = rbind(auc.sl.list,roc.sl.pert[[i]]$auc)
    auc.ssl.list = rbind(auc.ssl.list,roc.ssl.pert[[i]]$auc)
    auc.c.list = rbind(auc.c.list,roc.c.pert[[i]]$auc)
    roc.sl.list = rbind(roc.sl.list,roc.sl.pert[[i]]$roc)
    roc.ssl.list = rbind(roc.ssl.list,roc.ssl.pert[[i]]$roc)
    roc.c.list = rbind(roc.c.list,roc.c.pert[[i]]$roc)
  }
  
  # sl.tab = cbind( rbind(auc.sl, sd(auc.sl.list, na.rm=T)), 
  #                 rbind(roc.table.sl, apply(roc.sl.list[complete.cases(roc.sl.list),], 2, sd)))
  # rownames(sl.tab)=c('Pt.est','sd')
  # ssl.tab = cbind( rbind(auc.ssl, sd(auc.ssl.list, na.rm=T)), 
  #                  rbind(roc.table.ssl, apply(roc.ssl.list[complete.cases(roc.ssl.list),], 2, sd)))
  # rownames(ssl.tab)=c('Pt.est','sd')
  # true.tab = c(auc.true, roc.table.true)
  
  ## remove rows with NA in the bootstrap
  id1 = union(which(rowSums(is.na(auc.sl.list))==1), which(rowSums(is.na(auc.ssl.list))==1))
  id = union(id1, which(rowSums(is.na(auc.c.list))==1))
  auc.sl.list = auc.sl.list[-id,]
  auc.ssl.list = auc.ssl.list[-id,]
  auc.c.list = auc.c.list[-id,]
  
  se.auc.sl=sd(auc.sl.list, na.rm=T)
  se.roc.sl= tryCatch(apply(roc.sl.list[complete.cases(roc.sl.list),], 2, sd), error=function(e) NA)
  se.auc.ssl = sd(auc.ssl.list, na.rm=T)
  se.roc.ssl =  tryCatch(apply(roc.ssl.list[complete.cases(roc.ssl.list),], 2, sd), error=function(e) NA)
  se.auc.c = sd(auc.c.list, na.rm=T)
  se.roc.c =  tryCatch(apply(roc.c.list[complete.cases(roc.c.list),], 2, sd), error=function(e) NA)
  
  pt.auc.sl = mean(auc.sl.list, na.rm=T)
  pt.auc.ssl = mean(auc.ssl.list, na.rm=T)
  pt.auc.c = mean(auc.c.list, na.rm=T)
  
  return(list(auc.sl=pt.auc.sl, roc.table.sl=roc.table.sl, 
              auc.ssl=pt.auc.ssl, roc.table.ssl=roc.table.ssl,
              auc.c=pt.auc.c, roc.table.c=roc.table.c,
              auc.true =auc.true, roc.table.true =roc.table.true,
              se.auc.sl=se.auc.sl, se.roc.sl=se.roc.sl, 
              se.auc.ssl =se.auc.ssl, se.roc.ssl =se.roc.ssl,
              se.auc.c =se.auc.c, se.roc.c =se.roc.c))
}


analyze_result_sim = function(df.true, df.sl, df.ssl, df.comb, df.jess ){
  row = cbind(mean(df.true, na.rm = T), sd(df.true, na.rm = T),
              mean(df.sl, na.rm = T),   sd(df.sl, na.rm = T),
              mean(df.ssl, na.rm = T),  sd(df.ssl, na.rm = T),
              mean(df.comb, na.rm = T), sd(df.comb, na.rm = T),
              mean(df.jess, na.rm = T), sd(df.jess, na.rm = T))
  return(row)
}


