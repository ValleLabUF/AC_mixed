#calculates part of the marginal loglikelihood after integrating our the z's
get.calc.mloglik=function(dist.mat.sel,phi,dat,n.tsegm,n.ac,n.grid,theta){
  prob=exp(-phi*dist.mat.sel)
  prob1=prob/rowSums(prob)
  prob.mult=theta%*%prob1
  dat*log(prob.mult)
}
#-------------------------------------------------------------
sample.ac=function(ac.ind,dat,theta,n.ac,n.grid,phi,dist.mat,n.tsegm,n.possib.ac){
  ac.ind.orig=ac.ind.old=ac.ind

  for (i in 1:n.ac){
    ac.ind.new=ac.ind.old
    ind=sample(1:n.possib.ac,size=1)
    ac.ind.new[i]=ind
    
    #get marginal loglikel
    tmp.old=get.calc.mloglik(dist.mat.sel=dist.mat[ac.ind.old,],phi=phi,dat=dat,
                             n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,theta=theta)
    tmp.new=get.calc.mloglik(dist.mat.sel=dist.mat[ac.ind.new,],phi=phi,dat=dat,
                             n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,theta=theta)
    pold=sum(tmp.old)
    pnew=sum(tmp.new)
    
    #accept or reject MH
    k=acceptMH(p0=pold,p1=pnew,x0=ac.ind.old[i],x1=ac.ind.new[i],BLOCK=F)
    ac.ind.old[i]=k$x
  }
  list(ac.ind=ac.ind.old,accept=ac.ind.orig!=ac.ind.old)
}
#-----------------------------------
sample.phi=function(ac.ind,dist.mat,n.grid,n.ac,phi,jump,dat,n.tsegm,theta){
  old=phi
  new=abs(rnorm(1,mean=old,sd=jump)) #reflection proposal around zero
  
  #get marginal loglikel
  dist.mat1=dist.mat[ac.ind,]
  tmp.old=get.calc.mloglik(dist.mat.sel=dist.mat1,phi=old,dat=dat,
                           n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,theta=theta)
  tmp.new=get.calc.mloglik(dist.mat.sel=dist.mat1,phi=new,dat=dat,
                           n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,theta=theta)
  pold=sum(tmp.old)
  pnew=sum(tmp.new)
  
  #accept or reject MH
  k=acceptMH(p0=pold,p1=pnew,x0=old,x1=new,BLOCK=F)  
  logl=ifelse(k$accept==1,pnew,pold)
  list(phi=k$x,accept=k$accept)
}
#-----------------------------------
sample.z=function(ac.ind,dist.mat,n.grid,n.ac,n.tsegm,dat,phi,theta){
  #get distance and distance-related probabilities
  dist1=dist.mat[ac.ind,]
  prob=exp(-phi*dist1)
  prob.dist=prob/rowSums(prob)
  
  #sample z
  # z=array(0,dim=c(n.tsegm,n.grid,n.ac))
  # for (i in 1:n.tsegm){
  #   for (j in 1:n.grid){
  #     if (dat[i,j]>0){
  #       tmp=theta[i,]*prob.dist[,j]
  #       tmp1=tmp/sum(tmp)
  #       z[i,j,]=rmultinom(1,size=dat[i,j],prob=tmp1)
  #     }
  #   }
  # }
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
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.01
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
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