gibbs.activity.center=function(dat,grid.coord,n.ac,ac.coord.init,gamma1,possib.ac){
  #basic setup
  n.tsegm=nrow(dat)
  n.grid=nrow(grid.coord)
  grid.coord=data.matrix(grid.coord)
  n.possib.ac=nrow(possib.ac)
  n=rowSums(dat)
  
  #initial values
  ac.ind=sample(n.grid,size=n.ac)
  z=array(NA,dim=c(n.tsegm,n.grid,n.ac))
  for (i in 1:n.tsegm){
    for (j in 1:n.grid){
      z[i,j,]=rmultinom(1,size=dat[i,j],prob=rep(1/n.ac,n.ac))      
    }
  }
  phi=0.0001 #distance decay parameter
  theta=matrix(1/n.ac,n.tsegm,n.ac)
  
  #matrices to store results
  store.coord=matrix(NA,ngibbs,n.ac*2)
  store.param=matrix(NA,ngibbs,1) #to store phi
  store.logl=matrix(NA,ngibbs,1)
  store.theta=matrix(NA,ngibbs,n.ac*n.tsegm)
  
  #MH stuff
  adaptMH=50
  jump1=list(phi=0.2,ac.ind=rep(1,n.ac))
  accept1=list(phi=0,ac.ind=rep(0,n.ac))

  #pre-calculate distances between each potential AC location (possib.ac) and each actual location in our data (grid.coord)
  dist.mat=GetDistance(AcCoord=data.matrix(possib.ac),GridCoord=data.matrix(grid.coord), 
                       Ngrid=nrow(grid.coord), Nac=nrow(possib.ac))

  #gibbs sampler
  for (i in 1:ngibbs){
    print(i)
    
    #sample AC
    tmp=sample.ac(ac.ind=ac.ind,dat=dat,theta=theta,n.ac=n.ac,n.grid=n.grid,phi=phi,
                  dist.mat=dist.mat,n.tsegm=n.tsegm,n.possib.ac=n.possib.ac)
    ac.ind=tmp$ac.ind
    accept1$ac.ind=accept1$ac.ind+tmp$accept
    # ac.ind=ac.ind.true
    
    #sample phi
    # tmp=sample.phi(ac.ind=ac.ind,dist.mat=dist.mat,n.grid=n.grid,
    #                n.ac=n.ac,phi=phi,jump=jump1$phi,dat=dat,n.tsegm=n.tsegm,theta=theta)
    # phi=tmp$phi
    # accept1$phi=accept1$phi+tmp$accept
    phi=SamplePhi(w=0.02,MaxIter=100,phi=phi,dat=dat,theta=theta,DistMatSel=dist.mat[ac.ind,])
    # phi=phi.true
    
    #sample z
    z=sample.z(ac.ind=ac.ind,dist.mat=dist.mat,n.grid=n.grid,n.ac=n.ac,n.tsegm=n.tsegm,
               dat=dat,phi=phi,theta=theta)
    # z=z.true.disagg
    
    #sample theta
    v=sample.v(z=z,n.ac=n.ac,gamma1=gamma1,n.tsegm=n.tsegm)
    theta=GetTheta(v=v,nac=n.ac,ntsegm=n.tsegm)
    
    #get loglikel
    logl=GetCalcMloglik(dist.mat.sel=dist.mat[ac.ind,],phi=phi,dat=dat,theta=theta)

    if (i<nburn & i%%adaptMH==0){
      #adapt MH
      # tmp=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=adaptMH)
      # jump1=tmp$jump1
      # accept1=tmp$accept1
      print(accept1$ac.ind/i)
      
      #re-order data from time to time according to theta (largest to smallest)
      theta.m=apply(theta,2,median)
      ordem=order(theta.m,decreasing=T)
      ac.ind=ac.ind[ordem]
      theta=theta[,ordem]
      z=z[,,ordem]
    }
    
    #store results
    store.coord[i,]=unlist(grid.coord[ac.ind,])
    store.param[i,]=phi
    store.logl[i,]=logl
    store.theta[i,]=theta
  }
  list(coord=store.coord,phi=store.param,logl=store.logl,theta=store.theta,z=z)  
}
