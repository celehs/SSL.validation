
######## presented semi-supervised method

ini.FUN = function(S, Y){
  dat = cbind(S,Y)
  id.v=which(is.na(Y)!=1)
  S.full = S
  S.sorted = unique(sort(c(S)))
  k = length(S.sorted) 
  n = nrow(dat) 
  
  Dmat = Imat = matrix(1,nrow=n,ncol=k)
  Imat[1:n,1:k]=1*(S==VTM(S.sorted,n))
  
  Y.fake = g.logit(cbind(1,S.full)%*%glm(Y[id.v]~S.full[id.v], family="binomial", maxit = 10000)$coefficients)
  Y.fake[is.na(Y.fake)] = mean(Y[id.v]) 
  alpha1 = matrix(t(Y.fake)%*%Imat/t(Y.fake)%*%Dmat,nrow=k);
  alpha0 = matrix(t(1-Y.fake)%*%Imat/t(1-Y.fake)%*%Dmat,nrow=k);
  t = cbind(alpha1, alpha0)
  t
}

Est.EM = function(S, Y){
  dat = cbind(S,Y)
  id.v=which(is.na(Y)!=1)
  dat.v=dat[id.v,]
  dat.nv=dat[-id.v,]
  Yi.v = Y[id.v]
  Si.v = S[id.v]
  Si.nv = S[-id.v]
  Si = c(Si.v, Si.nv)
  p.ini = mean(Yi.v)
  
  Yihat.v = Yi.v
  
  Si.sorted = unique(sort(c(S)))
  k = length(Si.sorted)
  nn = nrow(dat)
  
  Dmat = Imat = matrix(1,nrow=nn,ncol=k)
  Imat.nv=matrix(1,nrow=dim(dat.nv)[1],ncol=k)
  Imat[1:nn,1:k]=1*(Si==VTM(Si.sorted,nn))
  Imat.nv[1:dim(dat.nv)[1],1:k]=1*(Si.nv==VTM(Si.sorted,dim(dat.nv)[1]))
  
  t.ini = ini.FUN(S,Y)
  alpha1.ini = t.ini[,1]
  alpha0.ini = t.ini[,2]
  
  p = p.ini
  alpha1 = alpha1.ini
  alpha0 = alpha0.ini
  k.step = 0
  eps=1;conv=1
  
  while(eps > 5e-5){
    alpha1.ini = alpha1
    alpha0.ini = alpha0
    p.ini = p
    
    k.step = k.step+1
    
    Yihat.nv = p.ini*Imat.nv%*%alpha1.ini/(p.ini*Imat.nv%*%alpha1.ini+(1-p.ini)*Imat.nv%*%alpha0.ini)
    Yihat.nv[Imat.nv%*%alpha1.ini==0]=0
    Yihat.nv[Imat.nv%*%alpha0.ini==0]=1
    
    Yihat = c(Yi.v, Yihat.nv)
    p = mean(Yihat)
    alpha1 = matrix(t(Yihat)%*%Imat/t(Yihat)%*%Dmat,nrow=k);
    alpha0 = matrix(t(1-Yihat)%*%Imat/t(1-Yihat)%*%Dmat,nrow=k);
    
    eps = sum(abs(c(p-p.ini, alpha1-alpha1.ini,alpha0-alpha0.ini)))
    if(is.na(eps)){conv=0;break}
    if(k.step>4000){conv=0;break}
  }
  print(k.step)
  if(conv==1){res=list(p=p, alpha1=alpha1, alpha0=alpha0, Si.sorted=Si.sorted)}
  if(conv==0){res=NA}
  res
}

ROC.est=function(alp1, alp0){
  sens.0=1-cumsum(alp1)
  omspec.0=1-cumsum(alp0)
  sens=rev(cumsum(alp1))
  omspec=rev(cumsum(alp0))
  sens=c(sens,0); omspec=c(omspec,0)
  height = (sens[-1]+sens[-length(sens)])/2
  width = -diff(omspec) 
  AUC=1-sum(height*width)
  list(AUC=AUC, fpr=omspec.0, tpr=sens.0)
}

#' Transform the data using quantiles and apply EM algorithm, return the auc and full ROC table (semi-supervised learning)
#' @param S surrogate S
#' @param Y labels containing NA
#' @return list with auc and roc table
#' @export
roc.semi.superv=function(S, Y){
  alpha=0.9
  h = floor(length(unique(S))^(alpha))+1
  dat.bin = dat.new.FUN(S,Y,h)
  dat.quant = dat.quantile.transf(dat.bin[,"S.new"],dat.bin[,"Y"])
  S = dat.quant[,"S.new"]
  Y = dat.quant[,"Y"]
  
  par.est=Est.EM(S, Y)
  alp1=par.est$alpha1
  alp0=par.est$alpha0
  p1=par.est$p; p0=1-p1
  cuts.t=par.est$Si.sorted
  cuts=quantile(dat.bin[,"S.new"],c(1:length(cuts.t))/length(cuts.t))
  
  junk=ROC.est(alp1, alp0) ## revise this so it returns AUC, fpr and tpr
  auc=junk$AUC
  tpr=junk$tpr
  fpr=junk$fpr
  ppv = tpr*p1/(tpr*p1+fpr*p0)
  npv = (1-fpr)*p0/((1-fpr)*p0+(1-tpr)*p1)
  p.pos=unlist(lapply(1:length(cuts), function(ll) mean(S>=cuts[ll])))
  junk =cbind("cut"= cuts,"p.pos"=p.pos, "fpr"=round(fpr,2),"tpr"=tpr,"ppv"=ppv,"npv"=npv)
  fpr0=seq(.01,.99,by=.01)
  id.print=unlist(lapply(fpr0, function(x) which.min(abs(x-fpr))[1]))
  out=list(auc=auc,roc=junk[id.print,],alpha1=alp1,alpha0=alp0)
  out
}

#' Apply a classic supervised method to return the AUC and full ROC table 
#' @param S surrogate of the small training set (all labeled)
#' @param Y labels of the small labeled set (no NA)
#' @return list with auc and roc table
#' @export
roc.superv=function(S,Y)
{
  yyi = S
  Di = Y
  yy0=0.5;fpr0=seq(0.01,0.99,0.01);wgti=NULL;yes.smooth=F
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
  if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
  mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
  for(k in 1:pp)
  {
    yy = yy0; 
    if(!is.null(fpr0)){
      tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
      TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y; 
      TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR); 
      yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))           
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }else{
      TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }
    out.yy = cbind(out.yy, yy)
    out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
    out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
    PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
    out.AUC <- c(out.AUC, AUC)
  }
  out.roc <- matrix(c(out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV), ncol=6)[-1,]
  colnames(out.roc) <- c("cut","p.pos","fpr","tpr","ppv","npv")
  out = list(auc=out.AUC,roc=out.roc)
  out
}


####### useful functions

g.logit = function(xx){exp(xx)/(exp(xx)+1)}

auc.FUN=function(TPR, FPR){
  sens=c(sort(TPR,decreasing=T),0); omspec=c(sort(FPR, decreasing=T),0)
  height = (sens[-1]+sens[-length(sens)])/2
  width = -diff(omspec) 
  AUC = sum(height*width)
  AUC
}

VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

S.FUN <- function(yy,Yi,Di,yes.smooth=F)
{
  if(yes.smooth){
    Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
    c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  }else{
    return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  }
  ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)
  ## sum_i I(yy FUN Yi)Vi
  # Vi weight
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
{
  yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth) 
  return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}
