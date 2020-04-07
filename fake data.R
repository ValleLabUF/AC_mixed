rm(list=ls(all=TRUE))
set.seed(111)

#basic setup
n.tsegm=400
n.ac=8
n.grid=100
n=floor(runif(n.tsegm,min=10,max=900))

#spatial coordinates
grid.coord=data.frame(x=runif(n.grid,min=0,max=100),
                      y=runif(n.grid,min=0,max=100))
ac.ind.true=ind=sample(nrow(grid.coord),size=n.ac)
ac.coord.true=ac.coord=grid.coord[ind,]

rangox=range(c(ac.coord$x,grid.coord$x))
rangoy=range(c(ac.coord$y,grid.coord$y))
plot(ac.coord$x,ac.coord$y,pch=19,xlim=rangox,ylim=rangoy)
points(grid.coord$x,grid.coord$y,col='red')

#generate theta
theta=matrix(NA,n.tsegm,n.ac)
for (i in 1:n.tsegm){
  if (i< n.tsegm/2){ #pure thetas
    theta[i,]=0
    ind=sample(1:n.ac,size=1)
    theta[i,ind]=1
  }
  if (i>= n.tsegm/2){
    tmp=runif(n.ac)
    theta[i,]=tmp/sum(tmp)
  }
}
image(theta)
theta.true=theta

#distance decay parameter
phi.true=phi=0.06

#calculate probabilities associated with each n.ac
probs=matrix(NA,n.ac,n.grid)
for (i in 1:n.ac){
  x2=(ac.coord$x[i]-grid.coord$x)^2
  y2=(ac.coord$y[i]-grid.coord$y)^2
  dist=sqrt(x2+y2)
  tmp=exp(-phi*dist)
  probs[i,]=tmp/sum(tmp)
}

#cluster membership and data
z=matrix(NA,n.tsegm,n.ac)
z.disagg=array(NA,dim=c(n.tsegm,n.grid,n.ac))
y=matrix(0,n.tsegm,n.grid)
for (i in 1:n.tsegm){
  z[i,]=rmultinom(1,size=n[i],prob=theta[i,])
  for (j in 1:n.ac){
    z.disagg[i,,j]=rmultinom(1,size=z[i,j],prob=probs[j,])
    y[i,]=y[i,]+z.disagg[i,,j]
  }
}
z.true.disagg=z.disagg
image(y)

setwd('U:\\GIT_models\\AC_mixed')
write.csv(y,'fake data.csv',row.names=F)
write.csv(grid.coord,'fake data grid.csv',row.names=F)