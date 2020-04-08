plot(res$logl,type='l')

plot(res$phi,type='l')
abline(h=phi.true,col='red')

ntot=n.tsegm*n.ac*n.grid
z.estim=apply(res$z,c(1,3),sum)
z.true=apply(z.true.disagg,c(1,3),sum)

#find ordem
ordem=rep(NA,ncol(z.true))
options(warn=0)
for (i in 1:ncol(z.true)){
  tmp=rep(NA,ncol(z.estim))
  for (j in 1:ncol(z.estim)){
    tmp[j]=cor(cbind(z.estim[,j],z.true[,i]))[1,2]
  }
  ind=which(!is.na(tmp) & tmp==max(tmp,na.rm=T))
  print(max(tmp))
  ordem[i]=ind
}
length(unique(ordem))
plot(z.estim[,ordem],z.true)

#--------------------
#look at theta
n.ac=20
tmp=matrix(res$theta[ngibbs-1,],n.tsegm,n.ac)
boxplot(tmp)
theta.estim=tmp[,ordem]
rango=range(c(theta.true,theta.estim))
plot(theta.true,theta.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')
#--------------------
ac.coord=matrix(res$coord[ngibbs-1,],n.ac,2)
ac.coord[ordem,]
ac.coord.true

rango=range(c(ac.coord.true),ac.coord[ordem,])
plot(ac.coord.true$x,ac.coord[ordem,1],xlim=rango,ylim=rango)
lines(rango,rango,col='red')

plot(ac.coord.true$y,ac.coord[ordem,2],xlim=rango,ylim=rango)
lines(rango,rango,col='red')

#--------------------
#look at spatial distribution
plot(possib.ac$x,possib.ac$y)
# points(ac.coord[,1],ac.coord[,2],col='blue',cex=2)
points(ac.coord[ordem,1],ac.coord[ordem,2],col='purple',cex=1.5)
points(ac.coord.true$x,ac.coord.true$y,col='red',cex=0.8,pch=19)