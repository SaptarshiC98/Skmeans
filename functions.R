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
s.kmeans=function(X,M,itermax=30){
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
  list1=list(label,M)
  names(list1)=c('label','centroids')
  return(list1)
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
