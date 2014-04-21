
combinesearch = function(X, y, M,flag=200,print=F,inipoint,Tem=20,annealing=F){
  time = proc.time()
  p = dim(X)[[2]]
  n = dim(X)[[1]]
  first = 1
  kset = c(0)
  varlist = list(1:p)
  betalist = list(rep(0,p))
  rss = 100 * sum(y^2)
  bestset = c()
  rssset = c(rss)
  iter = 100
  statu = 0
  stage = 0
  ns = n
  while(length(kset)>0){
    tempvar = varlist[[1]]
    varlist = varlist[-1]
    tempk = kset[1]
    kset = kset[-1]
    if(first){
      first = 0
      order = order(inipoint,decreasing=T) 
      tempvar = tempvar[order]
    }
    if(length(tempvar)<ns){
      stage = stage + 1
      Tem = Tem * 0.8
    }
    ns = length(tempvar)
    if(tempk >= M){
      next
    }
    if(ns<M){
      break
    }
    if(ns<n){
      sumbeta = gliner(X[,tempvar],y)
      sumrss = RSS(X[,tempvar],y,sumbeta)
      if(sumrss>rss){
        next
      }
    }
    flag = flag - 1 
    if(flag==0){
      break
    }
    hatbeta = gliner(X[,tempvar[1:M]],y)
    temprss = RSS(X[,tempvar[1:M]],y,hatbeta)
    
    if(temprss < rss){
      rss = temprss
      bestset = tempvar[1:M]
    }
    if(print==T){
      cat(flag,temprss,rss,stage,"\n")
    }
    if((ns>n)&(length(kset)>2*flag)){
      statu = 1
    }
    if(ns < n){
      statu = 0
    }
    if(statu==1){
      next
    }
    if(ns==M){
      break
    }
    if(tempk>(min(M,ns-2))){
      next
    }
    kset = c(kset, tempk:(min(M,ns-2)))
    for(i in (tempk+1):((min(M,ns-2))+1)){
      if(sum(is.na(tempvar[-i]))){
        browser()
      }
      varlist = c(varlist, list(tempvar[-i]))
    }
  }
  return(list(hatset=bestset,rss=rss,time=time))
}

combinetest = function(X, y, M,flag=200,print=F,inipoint,Tem=20,annealing=F){
  time = proc.time()
  p = dim(X)[[2]]
  n = dim(X)[[1]]
  first = 1
  kset = c(0)
  varlist = list(1:p)
  betalist = list(rep(0,p))
  rss = 100 * sum(y^2)
  bestset = c()
  rssset = c(rss)
  iter = 100
  statu = 0
  stage = 0
  ns = n
  while(length(kset)>0){
    tempvar = varlist[[1]]
    varlist = varlist[-1]
    tempk = kset[1]
    kset = kset[-1]
    inibeta = betalist[[1]]
    betalist = betalist[-1]
    inibeta = inibeta[tempvar]
    fatherrss = rssset[1]
    rssset = rssset[-1]
    if(length(tempvar)<ns){
      stage = stage + 1
      Tem = Tem * 0.8
    }
    ns = length(tempvar)
    if(tempk >= M){
      next
    }
    if(ns<M){
      break
    }
    if(ns<n){
      sumbeta = gliner(X[,tempvar],y)
      sumrss = RSS(X[,tempvar],y,sumbeta)
      if(sumrss>rss){
        next
      }
    }
    flag = flag - 1 
    if(flag==0){
      break
    }
    if(first){
      tempbeta = subset_EMtest(X[,tempvar], y, M, flag=iter,print=F,style="classic",ini="other",inipoint=inipoint)[[2]]
      first = 0
      iter = 1
    }  
    else{
      tempbeta = subset_EMtest(X[,tempvar], y, M, flag=iter,print=F,style="classic",ini="other",inipoint=inibeta)[[2]]
    }
    castnum = floor((ns - M - 1)/2)
    temporder = order(abs(tempbeta),decreasing=T)
    temporder[1:(M)] = temporder[(M):1]
    tmp = rep(0,p)
    tmp[tempvar] = tempbeta
    tempbeta = tmp
    tempvar = tempvar[temporder]
    tempvar = tempvar[1:(ns-castnum)]
    
    hatbeta = gliner(X[,tempvar[1:M]],y)
    temprss = RSS(X[,tempvar[1:M]],y,hatbeta)
    
    if(temprss < rss){
      rss = temprss
      bestset = tempvar[1:M]
    }
    if(annealing==T){
      if(cast(temprss,fatherrss,Tem)){
        next
      }
    }
    if(print==T){
      cat(flag,temprss,rss,stage,"\n")
    }
    if((ns>n)&(length(kset)>2*flag)){
      statu = 1
    }
    if(ns < n){
      statu = 0
    }
    if(statu==1){
      next
    }
    if(ns==M){
      break
    }
    if(tempk>(min(M,ns-2))){
      next
    }
    kset = c(kset, tempk:(min(M,ns-2)))
    for(i in (tempk+1):((min(M,ns-2))+1)){
      if(sum(is.na(tempvar[-i]))){
        browser()
      }
      varlist = c(varlist, list(tempvar[-i]))
    }
    for(i in (tempk+1):(min(M,ns-2)+1)){
      betalist = c(betalist, list(tempbeta))
    }
    rssset = c(rssset, rep(temprss,length(tempk:(min(M,ns-2)))))
  }
  return(list(hatset=bestset,rss=rss,time=time))
}


waste_combineradius = function(X, y, M,flag=200,print=F,inipoint,Tem=20,annealing=F,radiusceil=10){
  time = proc.time()
  p = dim(X)[[2]]
  n = dim(X)[[1]]
  first = 1
  kset = c(0)
  varlist = list(1:p)
  betalist = list(rep(0,p))
  rss = 100 * sum(y^2)
  bestset = c()
  rssset = c(rss)
  iter = 100
  statu = 0
  stage = 0
  ns = p
  while(length(kset)>0){
    tempvar = varlist[[1]]
    varlist = varlist[-1]
    tempk = kset[1]
    kset = kset[-1]
    inibeta = betalist[[1]]
    betalist = betalist[-1]
    inibeta = inibeta[tempvar]
    fatherrss = rssset[1]
    rssset = rssset[-1]
    if(length(tempvar)<ns){
      stage = stage + 1
      Tem = Tem * 0.4
    }
    radius = stage + tempk
    ns = length(tempvar)
    if(tempk >= M){
      next
    }
    if(ns<M){
      break
    }
    if(ns<n){
      sumbeta = gliner(X[,tempvar],y)
      sumrss = RSS(X[,tempvar],y,sumbeta)
      if(sumrss>rss){
        next
      }
    }
    flag = flag - 1 
    if(flag==0){
      break
    }
    if(radius<radiusceil){
      if(first){
        tempbeta = subset_EMtest(X[,tempvar], y, M, flag=iter,print=F,style="classic",ini="other",inipoint=inipoint)[[2]]
        first = 0
        iter=1
      }  
      else{
        tempbeta = subset_EMtest(X[,tempvar], y, M, flag=iter,print=F,style="classic",ini="other",inipoint=inibeta)[[2]]
      }
      castnum = floor((ns - M - 1)/4)
      temporder = order(abs(tempbeta),decreasing=T)
      temporder[1:(M)] = temporder[(M):1]
      tmp = rep(0,p)
      tmp[tempvar] = tempbeta
      tempbeta = tmp
      tempvar = tempvar[temporder]
      tempvar = tempvar[1:(ns-castnum)]
    }
    
    hatbeta = gliner(X[,tempvar[1:M]],y)
    temprss = RSS(X[,tempvar[1:M]],y,hatbeta)
    
    if(temprss < rss){
      rss = temprss
      bestset = tempvar[1:M]
    }
    if(annealing==T){
      if(cast(temprss,fatherrss,Tem)){
        next
      }
    }
    if(print==T){
      cat(flag,temprss,rss,stage,"\n")
    }
    if((ns>n)&(length(kset)>2*flag)){
      statu = 1
    }
    if(ns < n){
      statu = 0
    }
    if(statu==1){
      next
    }
    if(ns==M){
      break
    }
    if(tempk>(min(M,ns-2))){
      next
    }
    kset = c(kset, tempk:(min(M,ns-2)))
    for(i in (tempk+1):((min(M,ns-2))+1)){
      if(sum(is.na(tempvar[-i]))){
        browser()
      }
      varlist = c(varlist, list(tempvar[-i]))
    }
    for(i in (tempk+1):(min(M,ns-2)+1)){
      betalist = c(betalist, list(tempbeta))
    }
    rssset = c(rssset, rep(temprss,length(tempk:(min(M,ns-2)))))
  }
  return(list(hatset=bestset,rss=rss,time=time))
}



