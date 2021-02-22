########## simulation

## Generates the example data
## supervised St,Yt
## supervised true Sv,Yv
## SSL Sv,Yv.na
## JESS Sv, St, Yt
gen_data = function(n, N, p){
  #b0 = c(-5, 2, 2, 0.5, rep(0, p-3))
  b0 = c(-5, 2, 2, 3, rep(0, p-3))
  X = mvrnorm(N, rep(0,p) , diag(p)) + rbinom(N, 3, 0.3)
  Y = rbinom(N, 1, expit(cbind(1,X)%*%b0))
  dat.all = cbind(Y,X)
  inds.lab = sample(1:N, n)
  Yt = Y[inds.lab]; St = X[inds.lab,3]; Sv = X[,3]
  Yv.na = Y; Yv.na[-inds.lab]<-NA
  return(list(Yt = Yt, Yv.na = Yv.na, St = St, Sv = Sv, Yv = Y))
}

gen_data_p = function(n, N, p, pp){
  #b0 = c(-pp, 2, 2, 0.5, rep(0, p-3))
  b0 = c(-pp, 2, 2, 2, rep(0, p-3))
  X = mvrnorm(N, rep(0,p) , diag(p)) + rbinom(N, 3, 0.3)
  Y = rbinom(N, 1, expit(cbind(1,X)%*%b0))
  dat.all = cbind(Y,X)
  inds.lab = sample(1:N, n)
  Yt = Y[inds.lab]; St = X[inds.lab,3]; Sv = X[,3]
  Yv.na = Y; Yv.na[-inds.lab]<-NA
  return(list(Yt = Yt, Yv.na = Yv.na, St = St, Sv = Sv, Yv = Y))
}

sim.dat.FUN = function(nv,N,pp){
  mu = c(0,0.5,1,0.4)
  sigma = 5*diag(length(mu))
  beta0=  c(-5,pp,0.8,0.7,0.3)
  n = length(mu) 
  u = mvrnorm(N,mu,sigma)
  x = round(log(1+exp(u))) # we set initial beta as vector of ones ADD intercept ONLY CBIND not in x directly and add directly in beta
  pp = g.logit(cbind(1,x)%*%beta0)   # pp vector of size N
  Y = rbinom(N, 1, pp)
  id.t = sample(1:N, nv)
  Yt = Y[id.t] 
  Yv.na = Y; Yv.na[-id.t]<-NA
  x.t = x[id.t,]
  dt = as.data.frame.matrix(cbind(Yt,x.t))
  model = glm(Yt~x.t, family=binomial, data=dt)
  beta.hat = model$coefficients  # should be close to beta0
  S = cbind(1,x)%*%beta.hat # add intercept
  S = g.logit(S)  # S as probabilities
  St = S[id.t]
  return(list(Yt = Yt, Yv.na = Yv.na, St = St, Sv = S, Yv = Y))
}



sim.FUN=function(nv,N){
  n=10000
  valid_size=25
  p.0=0.4
  lambda1.0=3
  lambda0.0=1
  
  Y=rbinom(n, 1, p.0)
  S1=rpois(n, lambda=lambda1.0)
  S0=rpois(n, lambda=lambda0.0)
  S=S1*Y+S0*(1-Y)
  dat0=cbind.data.frame(Y, S)
  id.v=sample(1:n, nv, replace=F)
  dat0
}

library(MASS); library(stats); library(methods); library(MASS)
expit = function(xx){exp(xx)/(1+exp(xx))}
logit = function(xx){log(x/(1-x))}






######## presented semi-supervised method

# computes initial values in EM algorithm
# INPUT
# S : surrogate S
# Y : labels containing NA
# OUTPUT
# probability estimates alpha1, alpha0 
ini.FUN.pert = function(S, Y, Wt){
  id.v=which(is.na(Y)!=1)
  S.sorted = unique(sort(c(S)))
  k = length(S.sorted) 
  n = nrow(dat) 
  
  Dmat = Imat = matrix(1,nrow=n,ncol=k)
  Imat[1:n,1:k]=1*(S==VTM(S.sorted,n))
  
  #Y.fake = g.logit(cbind(1,S)%*%glm((Y[id.v]*Wt)~S[id.v], family="binomial", maxit = 10000)$coefficients)
  Y.fake = g.logit(cbind(1,S)%*%glm(Y[id.v]~S[id.v], family="binomial", maxit = 10000)$coefficients)
  # Y.fake[is.na(Y.fake)] = mean(Y[id.v]*Wt)/mean(Wt)
  Y.fake[is.na(Y.fake)] = mean(Y[id.v])
  alpha1 = matrix(t(Y.fake)%*%Imat/t(Y.fake)%*%Dmat,nrow=k);
  alpha0 = matrix(t(1-Y.fake)%*%Imat/t(1-Y.fake)%*%Dmat,nrow=k);
  p = mean(Y.fake)
  t = cbind(alpha1, alpha0, p)
  t
}

# EM algorithm
# INPUT
# S : surrogate S
# Y : labels containing NA
# OUTPUT
# estimate p, alpha1, alpha0
Est.EM.pert = function(S, Y, Wt){
  dat = cbind(S,Y)
  id.v=which(is.na(Y)!=1)
  dat.v=dat[id.v,]
  dat.nv=dat[-id.v,]
  Yi.v = Y[id.v]
  Si.v = S[id.v]
  Si.nv = S[-id.v]
  Si = c(Si.v, Si.nv)

  #Yihat.v = Yi.v*Wt
  
  Si.sorted = unique(sort(c(S)))
  k = length(Si.sorted)
  nn = nrow(dat)
  
  Dmat = Imat = matrix(1,nrow=nn,ncol=k)
  Imat.nv=matrix(1,nrow=dim(dat.nv)[1],ncol=k)
  Imat[1:nn,1:k]=1*(Si==VTM(Si.sorted,nn))
  Imat.nv[1:dim(dat.nv)[1],1:k]=1*(Si.nv==VTM(Si.sorted,dim(dat.nv)[1]))
  
  t.ini = ini.FUN.pert(S, Y, Wt)
  alpha1 = t.ini[,1]
  alpha0 = t.ini[,2]
  p = mean(Yi.v*Wt)/mean(Wt)

  k.step = 0
  eps=1;conv=1
  
  while(eps > 5e-5){
    alpha1.ini = alpha1
    alpha0.ini = alpha0
    p.ini = p
    
    k.step = k.step+1
    
    Yihat.nv = (p.ini*Imat.nv)%*%alpha1.ini/((p.ini*Imat.nv)%*%alpha1.ini+((1-p.ini)*Imat.nv)%*%alpha0.ini)
    Yihat.nv[Imat.nv%*%alpha1.ini==0]=0
    Yihat.nv[Imat.nv%*%alpha0.ini==0]=1
    
    Yihat = c(Yi.v, Yihat.nv)
    Wt.all = c(Wt, rep(1,length(Yihat.nv)))
    p = mean(Yihat*Wt.all)/mean(Wt.all)
    alpha1 = matrix(t(Yihat*Wt.all)%*%Imat/t(Yihat*Wt.all)%*%Dmat,nrow=k);
    alpha0 = matrix(t((1-Yihat)*Wt.all)%*%Imat/t((1-Yihat)*Wt.all)%*%Dmat,nrow=k);
    
    eps = sum(abs(c(p-p.ini, alpha1-alpha1.ini,alpha0-alpha0.ini)))
    if(is.na(eps)){conv=0;break}
    if(k.step>4000){conv=0;break}
  }
  #print(k.step)
  if(conv==1){res=list(p=p, alpha1=alpha1, alpha0=alpha0, Si.sorted=Si.sorted)}
  if(conv==0){res=NA}
  res
}

# Compute AUC, FPR, TPR using alpha1, alpha0
# ROC.est=function(alp1, alp0){
#   sens.0=1-cumsum(alp1)
#   omspec.0=1-cumsum(alp0)
#   sens=rev(cumsum(alp1))
#   omspec=rev(cumsum(alp0))
#   sens=c(sens,0); omspec=c(omspec,0)
#   height = (sens[-1]+sens[-length(sens)])/2
#   width = -diff(omspec) 
#   AUC=1-sum(height*width)
#   list(AUC=AUC, fpr=omspec.0, tpr=sens.0)
# }


# Transform data and apply EM algorithm, return the auc and full ROC table
# INPUT
# S : surrogate S
# Y : labels containing NA
# OUTPUT
# list with auc and roc table

roc.semi.superv.pert=function(S, Y, Wt){
  S.ini<-S
  Y.ini<-Y
  
  alpha=0.9
  h = floor(length(unique(S))^(alpha))+1
  dat.bin = dat.new.FUN(S,Y,h)
  dat.quant = dat.quantile.transf(dat.bin[,"S.new"], dat.bin[,"Y"])
  S = dat.quant[,"S.new"]
  Y = dat.quant[,"Y"]
  
  par.est=Est.EM.pert(S, Y, Wt)
  alp1=par.est$alpha1
  alp0=par.est$alpha0
  p1=par.est$p; p0=1-p1
  cuts=par.est$Si.sorted
  
  junk=ROC.est(alp1, alp0) ## revise this so it returns AUC, fpr and tpr
  auc=junk$AUC
  tpr=junk$tpr
  fpr=junk$fpr
  ppv = tpr*p1/(tpr*p1+fpr*p0)
  npv = (1-fpr)*p0/((1-fpr)*p0+(1-tpr)*p1)
  
  fscore = 2*ppv*tpr/(ppv+tpr)
  p.pos=unlist(lapply(1:length(cuts), function(ll) mean(S>=cuts[ll])))
  cuts=quantile(S.ini,1-p.pos)
  
  junk =cbind("cut"= cuts,"p.pos"=p.pos,"fpr"=round(fpr,2),"tpr"=tpr,"ppv"=ppv,"npv"=npv,"F.score"=fscore)
  fpr0=seq(.01,.99,by=.01)
  id.print=unlist(lapply(fpr0, function(x) which.min(abs(x-fpr))[1]))
  #out=list(auc=auc,roc=junk[id.print,],alpha1=alp1,alpha0=alp0)
  out=list(auc=auc,roc=junk[id.print,])
  out
}


# cuts.t=par.est$Si.sorted
# p.pos.t=unlist(lapply(1:length(cuts.t), function(ll) mean(S>=cuts.t[ll])))
# cuts=1-quantile(S.ini,p.pos.t)
# p.pos=unlist(lapply(1:length(cuts), function(ll) mean(S.ini>=cuts[ll])))

#cuts=par.est$Si.sorted
#cuts=1-quantile(dat.bin[,"S.new"],p.pos)
#cuts=1-quantile(S.ini,p.pos)
## S, S.new
# first get CUT_new (for example 0.03), and related p.pos
# find use p.pos to find CUT from S, so that P(S>CUT)=p.pos
# need to validate using simulation

########
get_beta_params <- function(mu, var){
  alpha <- ((1 - mu)/var - 1/mu) * mu^2
  beta <- alpha*(1/mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}



######## classic supervised method

# 
# roc.superv=function(S,Y,wgti=NULL)
# {
#   yyi = S
#   Di = Y
#   yy0=0.5;fpr0=seq(0.01,0.99,0.01);yes.smooth=F
#   out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- out.fscore <- NULL
#   if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
#   mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
#   for(k in 1:pp)
#   {
#     yy = yy0; 
#     if(!is.null(fpr0)){
#       tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
#       fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
#       TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y; 
#       TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR); 
#       yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))           
#       FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
#     }else{
#       TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
#       FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
#     }
#     out.yy = cbind(out.yy, yy)
#     out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
#     out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
#     PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
#     fscore <- 2*PPV*TPR/(PPV+TPR)
#     out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
#     out.fscore <- cbind(out.fscore,fscore)
#     AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
#     out.AUC <- c(out.AUC, AUC)
#   }
#   out.roc <- matrix(c(out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV,out.fscore), ncol=7)[-1,]
#   colnames(out.roc) <- c("cut","p.pos","fpr","tpr","ppv","npv","F.score")
#   out = list(auc=out.AUC,roc=out.roc)
#   out
# }
# 
# 
# ####### useful functions
# 
# g.logit = function(xx){exp(xx)/(exp(xx)+1)}
# 
# auc.FUN=function(TPR, FPR){
#   sens=c(sort(TPR,decreasing=T),0); omspec=c(sort(FPR, decreasing=T),0)
#   height = (sens[-1]+sens[-length(sens)])/2
#   width = -diff(omspec) 
#   AUC = sum(height*width)
#   AUC
# }
# 
# VTM<-function(vc, dm){
#   matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
# }
# 
# S.FUN <- function(yy,Yi,Di,yes.smooth=F)
# {
#   if(yes.smooth){
#     Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
#     c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
#   }else{
#     return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
#   }
#   ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
# }
# 
# sum.I <- function(yy,FUN,Yi,Vi=NULL)
#   ## sum_i I(yy FUN Yi)Vi
#   # Vi weight
# {
#   if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
#   # for each distinct ordered failure time t[j], number of Xi < t[j]
#   pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
#   if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
#   if (!is.null(Vi)) {
#     ## if FUN contains '=', tmpind is the order of decending
#     if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
#     ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
#     Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
#     return(rbind(0,Vi)[pos+1,])
#   } else return(pos)
# }
# # convert <- function(fit) {
# #   rochat.auc = fit$rochat[1]
# #   rochat.values = matrix(fit$rochat[-1],ncol=6)
# #   colnames(rochat.values) = c("cutoff","pos.rate","FPR","TPR","PPV","NPV")
# #   roc.cv.auc = fit$roc.cv[1]
# #   roc.cv.values = matrix(fit$roc.cv[-1],ncol=6)
# #   colnames(roc.cv.values) = c("cutoff","pos.rate","FPR","TPR","PPV","NPV")
# #   betahat = fit$beta[1:(grep("\\<bini",names(fit$beta))[1]-1)]
# #   names(betahat) = gsub("\\<b\\.","",names(betahat))
# #   return(list(ROC.hat.auc=rochat.auc, ROC.hat.values=rochat.values, ROC.cv.auc=roc.cv.auc,ROC.cv.values=roc.cv.values, beta=betahat, Y.hat=fit$Si))
# # }
# 
# Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
# {
#   yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth) 
#   return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
# }

coverage.calc.func = function(level=0.95, pt.est, se, gamma.true){
  ##
  ## pt.est and se are vectors of length=nsim
  ##
  ## level is confidence level, e.g. 0.95
  ci.band = qnorm(1-(1-level)/2)*se
  ci.left = pt.est - ci.band
  ci.right = pt.est + ci.band
  cover.logic = (ci.left <= gamma.true & ci.right >= gamma.true)
  cover.perc = sum(cover.logic, na.rm=TRUE)/length(cover.logic[is.na(cover.logic)==FALSE])
  lower.logic = ci.right < gamma.true
  lower.perc = sum(lower.logic, na.rm=TRUE)/length(lower.logic[is.na(lower.logic)==FALSE])
  upper.logic = ci.left > gamma.true
  upper.perc = sum(upper.logic, na.rm=TRUE)/length(upper.logic[is.na(upper.logic)==FALSE])
  ## return results
  return(list(cover.perc=cover.perc, lower.perc=lower.perc, upper.perc=upper.perc))
}





