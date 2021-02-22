
# Transform data and apply the combined method, return the auc and full ROC table
# INPUT
# S : surrogate S
# Y : labels containing NA
# OUTPUT
# list with auc and roc table
roc.SSL.combine=function(Sv, Yv.na, Wt = NULL, fpr = seq(0.01, 0.99, by = 0.01), percentile = F, bw = NULL){
  Y=Yv.na
  S=Sv
  
  S.ini<-S
  Y.ini<-Y
  n.v = length(Sv); n.t = length(Yt); if(is.null(Wt)){Wt = rep(1,n.t)};
  # whether or not to use percentile scale 
  #if(percentile == T){St = sum.I(St, '>=', Sv)/n.v; Sv = sum.I(Sv, '>=', Sv)/n.v}
  
  alpha=0.9
  h = floor(length(unique(S))^(alpha))+1
  dat.bin = dat.new.FUN(S,Y,h)
  dat.quant = dat.quantile.transf(dat.bin[,"S.new"],dat.bin[,"Y"])
  S = dat.quant[,"S.new"]
  Y = dat.quant[,"Y"]

  par.est=Est.EM(S, Y)
  alp1=par.est$alpha1
  alp0=par.est$alpha0
  p1.1=par.est$p; p0.1=1-p1.1
  cuts=par.est$Si.sorted
  
  junk=ROC.est(alp1, alp0) 
  #auc=junk$AUC
  tpr.1=junk$tpr
  fpr.1=junk$fpr
  p.pos=unlist(lapply(1:length(cuts), function(ll) mean(S>=cuts[ll])))
  cuts=quantile(S.ini,1-p.pos)
  nc = length(cuts)
  
  ### need to change St Yt Sv because of quantile data
  id.v = which(!is.na(Y))
  Yt = Y[id.v]
  St = S[id.v]
  Sv = S
  Yv.na = Y
  # smoothing for imputation
  if(is.null(bw)){bw = sd(St)/n.t^0.45}
  mhat = NP.REG(St, Yt, Sv,  bw, Wt = Wt) 
  p1.2=mean(mhat); p0.2 = 1-p1.2;
  # TPR, FPR, PPV, NPV at the observed values of the score
  #cuts = sort(Sv); 
  tpr.2=sum.I(cuts, "<=", Sv, mhat) /sum(mhat);
  fpr.2=sum.I(cuts, "<=", Sv, (1-mhat))/sum((1-mhat)); 

  ###### mean between the two methods
  tpr = (tpr.1+tpr.2)/2
  fpr = (fpr.1+fpr.2)/2
  p1 = (p1.1+p1.2)/2
  p0 = (p0.1+p0.2)/2
  
  ppv = tpr*p1/(tpr*p1+fpr*p0)
  npv = (1-fpr)*p0/((1-fpr)*p0+(1-tpr)*p1)
  fscore = 2*ppv*tpr/(ppv+tpr)
  auc = sum(tpr[-1]*(fpr[-nc]-fpr[-1]));
  
  #junk =cbind("cut"= cuts,"p.pos"=p.pos, "fpr"=round(fpr,2),"tpr"=tpr,"ppv"=ppv,"npv"=npv)
  junk =cbind("cut"= cuts,"p.pos"=p.pos,"fpr"=round(fpr,2),"tpr"=tpr,"ppv"=ppv,"npv"=npv,"F.score"=fscore)
  fpr0=seq(.01,.99,by=.01)
  id.print=unlist(lapply(fpr0, function(x) which.min(abs(x-fpr))[1]))
  #out=list(auc=auc,roc=junk[id.print,],alpha1=alp1,alpha0=alp0)
  out=list(auc=auc,roc=junk[id.print,])
  out
}





# 
# analysis_run = function(Sv, St, Yt, Wt = NULL, fpr = seq(0.01, 0.99, by = 0.01), percentile = T,
#                              pert = TRUE, prop = 1, b = 500, my.fpr=0.05){
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
#   roc.ssl.ob = ROC.FUN.SSL.NP(Sv, St, Yt, Wt, fpr = fpr, percentile = percentile)
#   roc.ssl = roc.ssl.ob$outpt[roc.ssl.ob$outpt[, 'FPR'] == my.fpr , ]
#   bw = roc.ssl.ob$bw
#   auc.ssl=roc.ssl.ob$auc
#   
#   ## run the resampling 
#   set.seed(1234)
#   
#   if(pert == TRUE){
#     
#     ## perturbation 
#     params = get_beta_params(0.25, (1/prop)/16)
#     inds.resamp = sapply(1:b, function(jj) 4*rbeta(n.t,params$alpha, params$beta))
#     
#     n = nrow(inds.resamp); bw.r = bw/(prop^0.45)
#     roc.ssl.pert = lapply(1:b, function(jj) 
#       tryCatch(ROC.FUN.SSL.NP(Sv, St, Yt, Wt*inds.resamp[,jj], fpr = fpr, percentile = percentile, bw = bw.r)$outpt), error=function(e) NA)
#     for (i in 1:b){roc.ssl.pert[[i]] = roc.ssl.pert[[i]][roc.ssl.pert[[i]][,'FPR']==my.fpr,]}
#     
#   }else{
#     ## bootstrap/subsampling
#     if(prop == 1){yes.replace = TRUE}else{yes.replace = FALSE}
#     inds.resamp = sapply(1:b, function(jj) sample(1:n.t, size = floor(prop*n.t), replace = yes.replace))
#     n = nrow(inds.resamp); bw.r = bw/(prop^0.45)
#     roc.ssl.pert = lapply(1:b, function(jj) {inds.tmp = inds.resamp[,jj]; 
#     tryCatch(ROC.FUN.SSL.NP(Sv, St[inds.tmp], Yt[inds.tmp], Wt[inds.tmp], fpr = fpr, percentile = percentile, bw = bw.r), error=function(e) NA)})
#     for (i in 1:b){roc.ssl.pert[[i]]$outpt = tryCatch(unlist(cutoff.choose(roc.ssl.pert[[i]]$outpt,"FPR",0.05)[2]), error=function(e) NA)}
#   }
#   
#   return(list(roc.ssl = roc.ssl.ob, roc.ssl.pert = roc.ssl.pert))
# }
# 
