library(igraph)
library(rgl)
library(MASS)
library(mlbench)
#DEFINE THE DELTA DISTANCE

delta=function(x,y){     ## x,y > 0 are vectors
  
  sqrt(  sum( log( (x+y)/(2*sqrt(x*y)))))
  
}



delta=function(x,y){     ## x,y > 0 are vectors
  
  sqrt( sum( log((x+y)/2)-0.5*(log(x)+log(y))))
  
}

# DEFINE WHAT THE CENTER IS
E=function(mu,Y){
  n=dim(Y)[1]
  s=0
  for(i in 1:n){
    s=s+(delta(Y[i,],mu))^2
  }
  return(s)
}

# KMEANS
s.kmeans=function(X,M,itermax){
  n=dim(X)[1]
  d=dim(X)[2]
  k=dim(M)[1]
  label=numeric(n)
  dist=numeric(k)
  
  for(l in 1: itermax){
    for(i in 1:n){
      for(j in 1:k){
        dist[j]=delta(X[i,],M[j,])
      }
      label[i]=which.min(dist)
    }
    
    
    for(i in 1:k){
      I=which(label==i)
      M[i,]=optim(colMeans(X[I,]),E,Y=X[I,])$par
    }
  }
  return(list(label,M))
}


euc.dist.sq=function(x1,x2){
  p=(x1-x2)^2
  return(sum(p))
}


k.means=function(X,M,tmax){
  
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  
  n=dim(X)[1]
  d=dim(X)[2]
  c=dim(M)[1]
  label=numeric(n)
  dist=numeric(c)
  t=0
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : c){
        dist[j]=euc.dist.sq(X[i,],M[j,])
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(i in 1:c){
      I=which(label==i)
      M[i,]=colMeans(X[I,])
    }
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M))
  
}

msq.euc.dist= function(x1, x2,a) (sum((abs(x1 - x2)) ^ a))^(1/a)

mwt.euc.dist.sq=function(x1,x2,a,w){
  p=(abs(x1-x2))^a
  p=w*p
  return(sum(p))
}

mvec.wt.euc.dist.sq=function(x1,x2,a,w){
  p=(abs(x1-x2))^a
  p=w*p
  return(p)
}

mE=function(mu,Y,b){
  n=dim(Y)[1]
  s=0
  for(i in 1:n){
    s=s+msq.euc.dist(Y[i,],mu,b)^b
  }
  return(s)
}


mwkmeans=function(X,M,a,beta,tmax){
  
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  
  n=dim(X)[1]
  d=dim(X)[2]
  c=dim(M)[1]
  weight=rep(1/d,d)
  label=numeric(n)
  dist=numeric(c)
  t=0
  D=numeric(d)
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : c){
        dist[j]=mwt.euc.dist.sq(X[i,],M[j,],a,weight^beta)
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(i in 1:c){
      I=which(label==i)
      M[i,]=optim(colMeans(X[I,]),mE,Y=X[I,],b=a)$par
    }
    
    #update weights
    for(j in 1:d){
      D[j]=0
    }
    for(i in 1:c){
      I=which(label==i)
      for(k in I){
        D=D+mvec.wt.euc.dist.sq(X[k,],M[i,],a,rep(1,d))
      }
    }
    
    for(i in 1:d){
      if(D[i]!=0){
        D[i]=1/D[i]
        D[i]=D[i]^(1/(beta-1))
      }
    }
    sum=sum(D)
    weight=D/sum
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M,weight))
  
}

sq.euc.dist= function(x1, x2) sum((x1 - x2) ^ 2)

wt.euc.dist.sq=function(x1,x2,w){
  p=(x1-x2)^2
  p=w*p
  return(sum(p))
}

vec.wt.euc.dist.sq=function(x1,x2,w){
  p=(x1-x2)^2
  p=w*p
  return(p)
}


wkmeans=function(X,M,beta,tmax){
  
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  
  n=dim(X)[1]
  d=dim(X)[2]
  c=dim(M)[1]
  weight=rep(1/d,d)
  label=numeric(n)
  dist=numeric(c)
  t=0
  D=numeric(d)
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : c){
        dist[j]=wt.euc.dist.sq(X[i,],M[j,],weight^beta)
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(i in 1:c){
      I=which(label==i)
      M[i,]=colMeans(X[I,])
    }
    
    #update weights
    for(j in 1:d){
      D[j]=0
    }
    for(i in 1:c){
      I=which(label==i)
      for(k in I){
        D=D+vec.wt.euc.dist.sq(X[k,],M[i,],rep(1,d))
      }
    }
    
    for(i in 1:d){
      if(D[i]!=0){
        D[i]=1/D[i]
        D[i]=D[i]^(1/(beta-1))
      }
    }
    sum=sum(D)
    weight=D/sum
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M,weight))
  
}



X=read.csv('Datasets/wall robot 2.csv',head=FALSE)
X=read.table('Datasets/mydata/sizea.txt')
data(iris)
X=iris
head(X)
class(X)
X=data.matrix(X)
toss=X[,4]
X=X[,-4]
n=dim(X)[1]
plot(X,col=toss)
num.iter=20
(numclus=max(toss))
sort(unique(toss))
indexing=matrix(rep(0,num.iter*8),ncol=8)
min(X)
X=X-min(X)+1
for(i in 1 : num.iter){
  #set.seed(100+i*77)
  sa=sample(n,numclus)
  M=X[sa,]
  #l1=kmeans(X,M,30)
  #l2=s.kmeans(X,M,30)
  l3=wkmeans(X,M,5,30)
  l4=mwkmeans(X,M,3.9,5,30)
  a1=l1[[1]]
  a2=l2[[1]]
  a3=l3[[1]]
  a4=l4[[1]]
  indexing[i,1]=compare(toss,a1,method='nmi')
  indexing[i,2]=compare(toss,a2,method='nmi')
  indexing[i,3]=compare(toss,a3,method='nmi')
  indexing[i,4]=compare(toss,a4,method='nmi')
  indexing[i,5]=compare(toss,a1,method='adjusted.rand')
  indexing[i,6]=compare(toss,a2,method='adjusted.rand')
  indexing[i,7]=compare(toss,a3,method='adjusted.rand')
  indexing[i,8]=compare(toss,a4,method='adjusted.rand')
  cat(i)
  cat('\n')
}
indexing
write.matrix(indexing, "results/density5_r.txt", sep="\t") 
colMeans(indexing)
setwd("C:/Users/User-PC/Desktop/ebooks/Sem 4/project/S divergence")
