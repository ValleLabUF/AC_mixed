#calculates part of the marginal loglikelihood after integrating our the z's
GetCalcMloglik=function(dist.mat.sel,phi,dat,theta){
  prob=exp(-phi*dist.mat.sel)
  prob1=prob/rowSums(prob)
  prob.mult=theta%*%prob1
  sum(dat*log(prob.mult))
}
#-------------------------------------------------------------
sample.ac=function(ac.ind,dat,theta,n.ac,n.grid,phi,dist.mat,n.tsegm,n.possib.ac){
  ac.ind.orig=ac.ind.old=ac.ind

  for (i in 1:n.ac){
    ac.ind.new=ac.ind.old
    ind=sample(1:n.possib.ac,size=1)
    ac.ind.new[i]=ind
    
    #get marginal loglikel
    pold=GetCalcMloglik(dist.mat.sel=dist.mat[ac.ind.old,],phi=phi,dat=dat,theta=theta)
    pnew=GetCalcMloglik(dist.mat.sel=dist.mat[ac.ind.new,],phi=phi,dat=dat,theta=theta)

    #accept or reject MH
    k=acceptMH(p0=pold,p1=pnew,x0=ac.ind.old[i],x1=ac.ind.new[i],BLOCK=F)
    ac.ind.old[i]=k$x
  }
  list(ac.ind=ac.ind.old,accept=ac.ind.orig!=ac.ind.old)
}
#-----------------------------------
# sample.phi=function(ac.ind,dist.mat,n.grid,n.ac,phi,jump,dat,n.tsegm,theta){
#   old=phi
#   new=abs(rnorm(1,mean=old,sd=jump)) #reflection proposal around zero
#   
#   #get marginal loglikel
#   dist.mat1=dist.mat[ac.ind,]
#   pold=GetCalcMloglik(dist.mat.sel=dist.mat1,phi=old,dat=dat,theta=theta)
#   pnew=GetCalcMloglik(dist.mat.sel=dist.mat1,phi=new,dat=dat,theta=theta)
# 
#   #accept or reject MH
#   k=acceptMH(p0=pold,p1=pnew,x0=old,x1=new,BLOCK=F)  
#   logl=ifelse(k$accept==1,pnew,pold)
#   list(phi=k$x,accept=k$accept)
# }
#-----------------------------------
sample.z=function(ac.ind,dist.mat,n.grid,n.ac,n.tsegm,dat,phi,theta){
  #get distance and distance-related probabilities
  dist1=dist.mat[ac.ind,]
  prob=exp(-phi*dist1)
  prob.dist=prob/rowSums(prob)
  
  #sample z
  tmp=SampleZ(ntsegm=n.tsegm, ngrid=n.grid, nac=n.ac,
              z=rep(0,n.tsegm*n.grid*n.ac),
              dat=dat, theta=theta, ProbDist=prob.dist)
  tmp$z
}
#-----------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually

  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1

  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)

  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]
  accept <- length(keep)

  list(x = x0, accept = accept)
}
#----------------------------
# print.adapt = function(accept1z,jump1z,accept.output){
#   accept1=accept1z; jump1=jump1z; 
#   
#   for (k in 1:length(accept1)){
#     z=accept1[[k]]/accept.output
#     print(names(accept1)[k])
#     print(mean(z)); print(mean(jump1[[k]]))
#   }
#   
#   for (k in 1:length(jump1)){
#     cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
#     jump1[[k]][cond] = jump1[[k]][cond]*2       
#     cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.01
#     jump1[[k]][cond] = jump1[[k]][cond]*0.5
#     accept1[[k]][]=0
#   }
#   
#   return(list(jump1=jump1,accept1=accept1))
# }
#----------------------------
sample.v=function(z,n.ac,gamma1,n.tsegm){
  ntk=apply(z,c(1,3),sum)

  seq1=n.ac:1
  tmp=t(apply(ntk[,seq1],1,cumsum))
  n.ge=tmp[,seq1]
  n.ge1=n.ge[,-1]
  
  v=rbeta(n.tsegm*(n.ac-1),ntk[,-n.ac]+1,n.ge1+gamma1)
  v1=matrix(v,n.tsegm,n.ac-1)
  cbind(v1,1)
}
#------------------------------------
#this function doubles the interval until we are outside the slice
DoublingPhi=function(yslice, w, phi, MaxIter, DistMatSel, dat, theta){

  ParamLo=phi-w*runif(1);
  ParamHi=ParamLo+w;
  ylo=GetCalcMloglik(dist.mat.sel=DistMatSel,phi=ParamLo,dat=dat,theta=theta);
  yhi=GetCalcMloglik(dist.mat.sel=DistMatSel,phi=ParamHi,dat=dat,theta=theta);
  
  oo=0;
  while((ylo>yslice) & (oo<MaxIter)){
    ParamLo=ParamLo-w;
    ylo=GetCalcMloglik(dist.mat.sel=DistMatSel,phi=ParamLo,dat=dat,theta=theta);
    oo=oo+1;
  }
  oo=0;
  while((yhi>yslice) & (oo<MaxIter)){
    ParamHi=ParamHi+w;
    yhi=GetCalcMloglik(dist.mat.sel=DistMatSel,phi=ParamHi,dat=dat,theta=theta)
    oo=oo+1;
  }
  c(ParamLo,ParamHi)
}
#--------------------------------------------------------
#this function shrinks the slice if samples are outside the slice. If sample is inside the slice, accept this sample
ShrinkPhi=function(rango1,yslice,dat,MaxIter,theta,DistMatSel) {
  yfim=-Inf;
  oo=0;
  diff1=Inf
  while ((yfim<yslice) & (diff1 > 0.00001) & (oo<MaxIter)){
    phi=runif(1,rango1[1],rango1[2])
    yfim=GetCalcMloglik(dist.mat.sel=DistMatSel,phi=phi,dat=dat,theta=theta);
    if (yfim<yslice){ #shrink the slice if phi falls outside
      DistLo=abs(rango1[1]-phi);
      DistHi=abs(rango1[2]-phi);
      if (DistLo<DistHi) rango1[1]=phi;
      if (DistLo>DistHi) rango1[2]=phi;
      diff1=rango1[2]-rango1[1];
    }
    oo=oo+1;
  }
  phi;
}
#---------------------------------------------------
#this function samples phi using a slice sampler 
SamplePhi=function(w,MaxIter,phi,dat,theta,DistMatSel) {

  #define upper bound
  upper1=GetCalcMloglik(dist.mat.sel=DistMatSel,phi=phi,dat=dat,theta=theta);
  yslice=upper1-rexp(1); #method suggest by Neal 2003 to sample uniformly vertically
      
  #define slice  
  rango1=DoublingPhi(yslice=yslice, w=w, phi=phi, MaxIter=MaxIter, 
                     DistMatSel=DistMatSel, dat=dat, theta=theta); #find range by doubling window
      
  #sample this particular parameter
  ShrinkPhi(rango1=rango1,yslice=yslice,dat=dat,MaxIter=MaxIter,
            theta=theta,DistMatSel=DistMatSel); #sample within the defined range (rango1)
}
