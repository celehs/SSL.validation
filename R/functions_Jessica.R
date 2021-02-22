## Semi-Supervised ROC Function
ROC.FUN.SSL.NP = function(Sv, St, Yt, Wt = NULL, fpr = seq(0.01, 0.99, by = 0.01), percentile = F, bw = NULL){
  
  ## St: score in the labeled set
  ## Yt: outcome
  ## Sv: score in the unlabeled set
  ## Wt: optional vector of weights
  ## fpr: desired fpr sequence for output
  ## percentile: indicator whether to use percentile scale
  
  n.v = length(Sv); n.t = length(Yt); if(is.null(Wt)){Wt = rep(1,n.t)};
  # whether or not to use percentile scale 
  if(percentile == T){St = sum.I(St, '>=', Sv)/n.v; Sv = sum.I(Sv, '>=', Sv)/n.v}
  # smoothing for imputation
  if(is.null(bw)){bw = sd(St)/n.t^0.45}
  mhat = NP.REG(St, Yt, Sv,  bw, Wt = Wt) 
  # prevalence
  mu1 = mean(mhat); mu0 = 1-mu1; 
  # TPR, FPR, PPV, NPV at the observed values of the score
  cuts = sort(Sv); nc = length(cuts); 
  TPR.c = sum.I(cuts, "<=", Sv, mhat) /sum(mhat);
  FPR.c = sum.I(cuts, "<=", Sv, (1-mhat))/sum((1-mhat)); 
  PPV.c = TPR.c*mu1/(TPR.c*mu1+FPR.c*mu0);
  NPV.c = (1-FPR.c)*mu0/((1-FPR.c)*mu0+(1-TPR.c)*mu1);
  fscore = 2*PPV.c*TPR.c/(PPV.c+TPR.c)
  # AUC
  auc = sum(TPR.c[-1]*(FPR.c[-nc]-FPR.c[-1]));
  # Interpolate to desired FPR sequence
  out = cbind("cut"= cuts, "p.pos"= mu1,"FPR"=FPR.c,"TPR"=TPR.c,"PPV"=PPV.c,"NPV"=NPV.c,"F.score"=fscore);
  out=sapply(1:ncol(out), function(kk){approx(out[,"FPR"], out[,kk], fpr, rule = 2, n = 100)$y}); 
  colnames(out) = c('cut', 'Prev', 'FPR', 'TPR', 'PPV', 'NPV', 'F.score')
  return(list(outpt = out, bw = bw, auc=auc))
}

SD.ROB = function(xx){
  ee = which(abs(xx-median(xx)) > 4*mad(xx)); 
  xx[ee] =NA; sd(xx, na.rm = T)}

SD.ROB.av = function(xx){(sd(xx) + mad(xx))/2}
## Non-parametric regression function
NP.REG = function(St, Yt, Sv, bw, Wt = NULL, kern.mat=NULL){
  nv = length(Sv); nt = length(St); if (is.null(Wt)){Wt = rep(1,nt)};
  if(is.null(kern.mat)){kern.mat = sapply(1:nt, function(kk) dnorm(Sv - rep(St[kk], nv), sd = bw))}
  nw.est =  kern.mat %*% (Yt*Wt) * (1/ (kern.mat%*%Wt)); 
  return(nw.est);
}


analysis_run.jess = function(Sv, St, Yt, Wt = NULL, fpr = seq(0.01, 0.99, by = 0.01), percentile = T,
                             pert = TRUE, prop = 1, b = 500, my.fpr=0.05){
  
  ## St: score in the labeled set
  ## Yt: outcome
  ## Sv: score in the unlabeled set
  ## Wt: optional vector of weights
  ## fpr: desired fpr sequence for output
  ## percentile: indicator whether to use percentile scale
  ## pert: indicator whether to use perturbation (TRUE) or bootstrap (FALSE)
  ## prop: proportion of data for resampling 
  ## b: number of resamples
  ## ex. for subsampling set pert = FALSE and prop = 0.6 to subsample 60% of data wo replacement
  
  n.t = length(Yt); if(is.null(Wt)){Wt = rep(1, n.t)}
  ## get the estimates for entire data
  roc.ssl.ob = ROC.FUN.SSL.NP(Sv, St, Yt, Wt, fpr = fpr, percentile = percentile)
  roc.ssl = roc.ssl.ob$outpt[roc.ssl.ob$outpt[, 'FPR'] == my.fpr , ]
  bw = roc.ssl.ob$bw
  auc.ssl=roc.ssl.ob$auc
  
  ## run the resampling 
  set.seed(1234)
  
  if(pert == TRUE){
    
    ## perturbation 
    params = get_beta_params(0.25, (1/prop)/16)
    inds.resamp = sapply(1:b, function(jj) 4*rbeta(n.t,params$alpha, params$beta))
    
    n = nrow(inds.resamp); bw.r = bw/(prop^0.45)
    roc.ssl.pert = lapply(1:b, function(jj) 
      tryCatch(ROC.FUN.SSL.NP(Sv, St, Yt, Wt*inds.resamp[,jj], fpr = fpr, percentile = percentile, bw = bw.r)$outpt), error=function(e) NA)
    for (i in 1:b){roc.ssl.pert[[i]] = roc.ssl.pert[[i]][roc.ssl.pert[[i]][,'FPR']==my.fpr,]}
    
  }else{
    ## bootstrap/subsampling
    if(prop == 1){yes.replace = TRUE}else{yes.replace = FALSE}
    inds.resamp = sapply(1:b, function(jj) sample(1:n.t, size = floor(prop*n.t), replace = yes.replace))
    n = nrow(inds.resamp); bw.r = bw/(prop^0.45)
    roc.ssl.pert = lapply(1:b, function(jj) {inds.tmp = inds.resamp[,jj]; 
    tryCatch(ROC.FUN.SSL.NP(Sv, St[inds.tmp], Yt[inds.tmp], Wt[inds.tmp], fpr = fpr, percentile = percentile, bw = bw.r), error=function(e) NA)})
    for (i in 1:b){roc.ssl.pert[[i]]$outpt = tryCatch(unlist(cutoff.choose(roc.ssl.pert[[i]]$outpt,"FPR",0.05)[2]), error=function(e) NA)}
  }

  return(list(roc.ssl = roc.ssl.ob, roc.ssl.pert = roc.ssl.pert))
}


analyze_results.jess <- function( roc.ssl, roc.ssl.pert, my.fpr = NULL, prop = 1, b=boot){
  
  if(prop != 1){fac = sqrt(1/(1 - prop))}else{fac = 1} # scaling for factor for finite population sampling 
  
  # m = ncol(roc.sl.pert); params.int = c(1, 3:7); 
  # 
  # sd.ssl <- sapply(params.int, function(kk) apply(roc.ssl.pert[, seq(kk, m, by = 7)], 1, sd))*fac
  # sd.rob.ssl <- sapply(params.int, function(kk) apply(roc.ssl.pert[, seq(kk, m, by = 7)], 1, SD.ROB))*fac
  # sd.rob.av.ssl <- sapply(params.int, function(kk) apply(roc.ssl.pert[, seq(kk, m, by = 7)], 1, SD.ROB.av))*fac
  # mean.ssl <- sapply(params.int, function(kk) apply(roc.ssl.pert[, seq(kk, m, by = 7)], 1, mean))
    
    auc.ssl = roc.ssl$auc
    roc.table.ssl = roc.ssl$outpt
    
    auc.ssl.list = NULL
    roc.ssl.list = NULL
    for (i in 1:b){
      auc.ssl.list = rbind(auc.ssl.list,roc.ssl.pert[[i]]$auc)
      roc.ssl.list = rbind(roc.ssl.list,roc.ssl.pert[[i]]$outpt)
    }
    se.auc.ssl = tryCatch(sd(auc.ssl.list, na.rm=T), error=function(e) NA)
    se.roc.ssl = tryCatch(apply(roc.ssl.list[complete.cases(roc.ssl.list),], 2, sd), error=function(e) NA)
    
    return(list( auc.ssl=auc.ssl, roc.table.ssl=roc.table.ssl,
                 se.auc.ssl =se.auc.ssl, se.roc.ssl =se.roc.ssl))
  
}

# ## Supervised ROC Function
# ROC.FUN.SL.NP = function(St, Yt, Wt = NULL, fpr = seq(0.01, 0.99, by = 0.01), percentile = T){
#   
#   ## St: score
#   ## Yt: outcome
#   ## Wt: optional vector of weights
#   ## fpr: desired fpr sequence for output
#   ## percentile: indicator whether to use percentile scale
#   
#   n.t = length(Yt); if(is.null(Wt)){Wt = rep(1,n.t)};
#   # whether or not to use percentile scale 
#   if(percentile == T){St = sum.I(St, '>=', St)/n.t}
#   # prevalence
#   mu1 = mean(Yt*Wt)/mean(Wt); mu0 = 1-mu1; 
#   # TPR, FPR, PPV, NPV at the observed values of the score
#   cuts = sort(St); nc = length(cuts); 
#   TPR.c = sum.I(cuts, "<=", St, Yt*Wt) /sum(Yt*Wt);
#   FPR.c = sum.I(cuts, "<=", St, (1-Yt)*Wt)/sum((1-Yt)*Wt); 
#   PPV.c = TPR.c*mu1/(TPR.c*mu1+FPR.c*mu0);
#   NPV.c = (1-FPR.c)*mu0/((1-FPR.c)*mu0+(1-TPR.c)*mu1);
#   # AUC
#   auc = sum(TPR.c[-1]*(FPR.c[-nc]-FPR.c[-1]));
#   # Interpolate to desired FPR sequence
#   out = cbind("cut"= cuts,"FPR"=FPR.c,"TPR"=TPR.c,"PPV"=PPV.c,"NPV"=NPV.c);
#   out=sapply(1:ncol(out), function(kk){approx(out[,"FPR"], out[,kk], fpr, rule = 2, n = 100)$y}); 
#   outpt = cbind(out, auc, mu1)
#   colnames(outpt) = c('cut', 'FPR', 'TPR', 'PPV', 'NPV', 'AUC', 'Prev')
#   return(outpt)
# }
## Fast summation function based on ranks
# ## Main function to run the analysis with pertrubation or bootstrap
# analysis_run.jess.original = function(Sv, St, Yt, Wt = NULL, fpr = seq(0.01, 0.99, by = 0.01), percentile = T,
#                         pert = TRUE, prop = 1, b = 500){
#   
#   ## St: score in the labeled set
#   ## Yt: outcome
#   ## Sv: score in the unlabeled set
#   ## Wt: optional vector of weights
#   ## fpr: desired fpr sequence for output
#   ## percentile: indicator whether to use percentile scale
#   ## pert: indicator whether to use perturbation (TRUE) or bootstrap (FALSE)
#   ## prop: proportion of data for resampling 
#   ## b: number of resamples
#   ## ex. for subsampling set pert = FALSE and prop = 0.6 to subsample 60% of data wo replacement
#   
#   n.t = length(Yt); if(is.null(Wt)){Wt = rep(1, n.t)}
#   ## get the estimates for entire data
#   roc.sl = ROC.FUN.SL.NP(St, Yt, Wt, fpr = fpr, percentile = percentile)
#   roc.ssl.ob = ROC.FUN.SSL.NP(Sv, St, Yt, Wt, fpr = fpr, percentile = percentile)
#   roc.ssl = roc.ssl.ob$outpt[my.fpr,]
#   bw = roc.ssl.ob$bw
#   
#   ## run the resampling 
#   set.seed(92047)
#   
#   if(pert == TRUE){
#     
#     ## perturbation 
#     params = get_beta_params(0.25, (1/prop)/16)
#     inds.resamp = sapply(1:b, function(jj) 4*rbeta(n.t,params$alpha, params$beta))
#     
#     n = nrow(inds.resamp); bw.r = bw/(prop^0.45)
#     roc.ssl.pert = lapply(1:b, function(jj) ROC.FUN.SSL.NP(Sv, St, Yt, Wt*inds.resamp[,jj], fpr = fpr, 
#                                                            percentile = percentile, bw = bw.r)$outpt)
#     roc.ssl.pert.rtn = do.call(cbind, roc.ssl.pert)
#     
#     roc.sl.pert = lapply(1:b, function(jj) ROC.FUN.SL.NP(St, Yt, Wt*inds.resamp[,jj], fpr = fpr,
#                                                          percentile = percentile))
#     roc.sl.pert.rtn = do.call(cbind, roc.sl.pert)
#     
#   }else{
#     
#     ## bootstrap/subsampling
#     if(prop == 1){yes.replace = TRUE}else{yes.replace = FALSE}
#     inds.resamp = sapply(1:b, function(jj) sample(1:n.t, size = floor(prop*n.t), replace = yes.replace))
#     
#     n = nrow(inds.resamp); bw.r = bw/(prop^0.45)
#     roc.ssl.pert = lapply(1:b, function(jj) {inds.tmp = inds.resamp[,jj]; 
#     ROC.FUN.SSL.NP(Sv, St[inds.tmp], Yt[inds.tmp], Wt[inds.tmp], fpr = fpr, percentile = percentile, bw = bw.r)$outpt})
#     roc.ssl.pert.rtn = do.call(cbind, roc.ssl.pert)
#     
#     roc.sl.pert = lapply(1:b, function(jj){inds.tmp = inds.resamp[,jj]; 
#     ROC.FUN.SL.NP(St[inds.tmp], Yt[inds.tmp], Wt[inds.tmp], fpr = fpr, percentile = percentile)})
#     roc.sl.pert.rtn = do.call(cbind, roc.sl.pert)
#   }
#   
#   return(list(roc.sl = roc.sl, roc.ssl = roc.ssl, roc.ssl.pert = roc.ssl.pert.rtn, roc.sl.pert =roc.sl.pert.rtn))
#   
# }
# 
# analyze_results.jess.original <- function(roc.sl, roc.ssl, roc.ssl.pert, roc.sl.pert, my.fpr = NULL, prop = 1){
#   
#   if(prop != 1){fac = sqrt(1/(1 - prop))}else{fac = 1} # scaling for factor for finite population sampling 
#   
#   m = ncol(roc.sl.pert); params.int = c(1, 3:7); 
#   sd.sl <- sapply(params.int, function(kk) apply(roc.sl.pert[, seq(kk, m, by = 7)], 1, sd))*fac
#   sd.rob.sl <- sapply(params.int, function(kk) apply(roc.sl.pert[, seq(kk, m, by = 7)], 1, SD.ROB))*fac
#   sd.rob.av.sl <- sapply(params.int, function(kk) apply(roc.sl.pert[, seq(kk, m, by = 7)], 1, SD.ROB.av))*fac
#   mean.sl <- sapply(params.int, function(kk) apply(roc.sl.pert[, seq(kk, m, by = 7)], 1, mean))
# 
#   sd.ssl <- sapply(params.int, function(kk) apply(roc.ssl.pert[, seq(kk, m, by = 7)], 1, sd))*fac
#   sd.rob.ssl <- sapply(params.int, function(kk) apply(roc.ssl.pert[, seq(kk, m, by = 7)], 1, SD.ROB))*fac
#   sd.rob.av.ssl <- sapply(params.int, function(kk) apply(roc.ssl.pert[, seq(kk, m, by = 7)], 1, SD.ROB.av))*fac
#   mean.ssl <- sapply(params.int, function(kk) apply(roc.ssl.pert[, seq(kk, m, by = 7)], 1, mean))
#   
#   if(is.null(my.fpr)){
#     return(list(roc.sl = roc.sl, roc.ssl = roc.ssl, sd.sl = sd.sl, sd.rob.sl = sd.rob.sl, sd.rob.av.sl = sd.rob.av.sl,
#                 sd.rob.ssl = sd.rob.ssl, sd.rob.av.ssl = sd.rob.av.ssl))
#   }else{
#     sl.tab = cbind(my.fpr, rbind(roc.sl[my.fpr, params.int], mean.sl[my.fpr,],
#                                  sd.sl[my.fpr,], sd.rob.sl[my.fpr, ], sd.rob.av.sl[my.fpr, ]))
#     ssl.tab = cbind(my.fpr, rbind(roc.ssl[my.fpr, params.int], mean.ssl[my.fpr,],
#                                   sd.ssl[my.fpr,], sd.rob.ssl[my.fpr, ], sd.rob.av.ssl[my.fpr, ]))
#     re.tab = cbind(my.fpr, (sl.tab[-c(1,2), 2:7]/ssl.tab[-c(1,2), 2:7])^2)
#     colnames(sl.tab) = colnames(ssl.tab) = colnames(re.tab) = c('FPR', 'cut', 'TPR', 'PPV', 'NPV', 'AUC', 'Prev')
#     rownames(sl.tab) = rownames(ssl.tab) = c('Pt Est', 'Mean Resamps.', 'SD', 'SD.ROB', '(SD+MAD)/2')
#     rownames(re.tab) = c('RE_SD', 'RE_SD.ROB', 'RE_(SD+MAD)/2')
#     return(list(sl.tab = sl.tab, ssl.tab = ssl.tab, re.tab = re.tab))
#   }
# }

