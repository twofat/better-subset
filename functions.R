require(MASS)
require(glmnet)
require(lars)
require(xtable)
require(ncvreg)
source("E://Variable Selective Methods//better subset//functions2.r")


tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}

generate = function(n,p,rho){
  sigma = matrix(,p,p)
  for(i in 1:p){
    for(j in 1:p){
      sigma[i,j] = rho^(abs(i-j))
    }
  }
  mvrnorm(n,rep(0,p),sigma)
}

gcvridge = function(y, X, num = 100){
  p = length(X[1,])
  n = length(X[,1])
  max = 10 * max(log(abs(X)))
  min = -max
  lambdaset = exp(seq(min,max,length.out=num))
  minV = 10^10;
  minlambda = 0;
  X_2 = t(X)%*%X
  for(lambda in lambdaset){
    A = X %*% solve(X_2 + lambda*diag(p)) %*% t(X)
    V = apply((diag(n) - A) %*% as.matrix(y),2, function(x) sum(x^2)) / ( sum(diag(diag(n)-A)) )^2   
    if(minV > V){
      minV = V
      minlambda = lambda
    }
  }  
  beta = as.vector(solve(t(X)%*%X + minlambda*diag(p)) %*% t(X) %*% as.matrix(y))
  return(list(beta,minlambda))
}

randomdesign = function(n,p,rho,truelength,sd){
  X = generate(n,p,rho)
  trueset = sample(1:p,truelength)
  truebeta = rep(0,p)
  truebeta[trueset] = rep(1,truelength)
  epsilon = rnorm(n,0,sd)
  y = X %*% as.matrix(truebeta) + epsilon 
  return(list(X=X,y=y,trueset=trueset))
}

S_M = function(x, M){
  if(is.matrix(x)){
    x = as.vector(x)
  }
  x = abs(x)
  result = order(x,decreasing=T)[1:M]  
  return(result)
}

delta_compute = function(X){
  p = dim(X)[2]
  if(n>p){
  temp = eigen(t(X)%*%X)
  }
  else{
    temp = eigen(X%*%t(X))
  }
  value = temp$value
  Q = temp$vectors
  a = max(temp$value)
  return(a)
}

gliner = function(X, y){
  y = as.matrix(y)
  temp1 = t(X)%*%y
  temp = ginv(t(X)%*%X)
  beta = temp %*% temp1
  return(beta)
}

beta_EM = function(X, inibeta, M, a, y){
  p = dim(X)[2]
  inibeta = as.matrix(inibeta)
  y = as.matrix(y)
  temp1 = t(X)%*%y/a
  temp2 = (diag(p)-t(X)%*%X/a)%*%inibeta
  temp = temp1 + temp2
  
  return(temp)
}
RSS = function(X, y, beta){
  y = as.matrix(y)
  beta = as.matrix(beta)
  return(sum((y-X%*%beta)^2))
}

subset_EM = function(X, y, M,print=F, flag=10,ini="zero",inipoint,style="classic"){
  time = proc.time()
  p = dim(X)[2]
  a = delta_compute(X)
  if(ini=="stepforward"){
    inipoint = stepforward(X,y,M)[[1]]
  }
  if(ini=="zero"){
    inipoint = rep(0,p)
    if(flag==0) flag=1
  }
  if(ini=="lar"){
    inipoint = coef(lars(X, y, "lar",intercept=F))[M+1,]
  }
  if(ini=="lasso"){
    lassofit = cv.glmnet(y=y,x=X,alpha=1,intercept=F)
    inipoint = coef(lassofit)[2:(p+1)]
    if(flag==0) flag=1
  }
  if(ini=="glm"){
    inipoint = is.numeric(tryCatch.W.E(gliner(X,y))$value)
    if(!is.numeric(inipoint)) inipoint = c(0,gliner(X[,-1],y))
    if(flag==0) flag=1
  }
  if(ini=="elastic"){
    if(flag==0) flag=1
    cv = c()
    betaen = c()
    temp = cv.glmnet(y=y,x=X,alpha=0.5,intercept=F)
    con = coef(temp)[1]
    for(lambda2 in 10^(-2:2)){
      X1temp = rbind(X, sqrt(lambda2)*diag(p))
      ytemp = c(y, rep(con,p))
      elasticfit = cv.glmnet(y=ytemp,x=X1temp,alpha=1,intercept=F)
      cv = c(cv, min(elasticfit$cvm))
      betaen = rbind(betaen, coef(elasticfit)[2:(p+1)])
    }    
    inipoint = betaen[which.min(cv),]        
  }
  if(ini=="scad"){
    if(flag==0) flag=1
    scadfit = cv.ncvreg(X,y,family="gaussian",penalty="SCAD")
    betascad = scadfit$fit$beta[,scadfit$min]
    inipoint = betascad[2:(p+1)]
  }
  inibeta = inipoint
  tempbeta = inibeta
  rss = RSS(X,y,inibeta)
  while(flag>0){
    flag = flag - 1
    if(style=="fast"){
      tempbeta = beta_EM(X=X,y=y,inibeta=inibeta,M=M,a=a)
      order = S_M(tempbeta, M)
      tempbeta = gliner(X[,order],y)
      inibeta = rep(0,p)
      inibeta[order] = tempbeta
    }
    if(style=="classic"){
      tempbeta = beta_EM(X=X,y=y,inibeta=inibeta,M=M,a=a)
      order = S_M(tempbeta, M)
      inibeta = rep(0,p)
      inibeta[order] = tempbeta[order]
    }
    oldrss = rss
    rss=sum((y-X%*%as.matrix(inibeta))^2)
    if(abs(rss-oldrss)<10^(-4)){
      break
    }
    if(print==T){
      print(rss)
    }
  }
  time = proc.time() - time
  hatset = sort((1:p)[inibeta!=0])
  return(list(inibeta,tempbeta,rss=rss,time=time,hatset=hatset))    
}

stepforward = function(X,y,M){
  time = proc.time()
  data = data.frame(X,y)
  null = lm(y~1,data)
  full = lm(y~.,data)
  result = stepAIC(null, steps=M,scope=list(lower=null,upper=full), direction="forward",trace=F)
  coef = result$coe[-1]
  set = sort(as.numeric(substring(names(coef),2)))
  iniset = set
  temp = lm(y~X[,iniset]+0)
  inibeta = rep(0,p)
  inibeta[iniset] = temp$coe
  rss = RSS(X[,iniset],y, temp$coe)
  time = proc.time() - time 
  hatset = iniset
  return(list(inibeta,rss=rss,time=time,hatset=hatset))
}

combine = function(X, y, M){
  p = dim(X)[[2]]
  n = dim(X)[[1]]
  kset = c(0)
  varlist = list(1:p)
  rss = 100 * sum(y^2)
  bestset = c()
  num = 0
  while(length(kset)>0){
    num = num + 1 
    print(num)
    tempvar = varlist[[1]]
    varlist = varlist[-1]
    tempk = kset[1]
    kset = kset[-1]
    ns = length(tempvar)
    tempvar = tempvar[sample(1:ns)]
    if(tempk >= M){
      next
    }
    if(ns<M){
      break
    }
    #    print(tempvar)
    print(tempvar[1:M])
    #    print(tempk)
    tempbeta = gliner(X[,tempvar[1:M]],y)
    temprss = RSS(X[,tempvar[1:M]],y,tempbeta)
    #    print(temprss)
    if(temprss < rss){
      rss = temprss
      bestset = tempvar[1:M]
    }
    kset = c(kset, tempk:(min(M,ns-2)))
    for(i in (tempk+1):((min(M,ns-2))+1)){
      varlist = c(varlist, list(tempvar[-i]))
    }
  }
  return(bestset)
}

cast = function(rssNew, rssMin, Tem){
  delta = rssNew - rssMin
  p = exp(-delta/Tem) 
  if(runif(1)>p){
    result = T
  } 
  else
    result = F
  if(delta<0)
    result = F
  return(result)
}

subset_EMtest = function(X, y, M,print=F, flag=10,ini="zero",inipoint,style="classic"){
  time = proc.time()
  p = dim(X)[2]
  a = delta_compute(X)
  rss = 100*sum(y^2)
  if(ini=="stepforward"){
    inipoint = stepforward(X,y,M)[[1]]
  }
  if(ini=="zero"){
    inipoint = rep(0,p)
  }
  inibeta = inipoint
  for(i in 1:flag){
    if(style=="fast"){
      tempbeta = beta_EM(X=X,y=y,inibeta=inibeta,M=M,a=a)
      order = S_M(tempbeta, M)
      tempbeta = gliner(X[,order],y)
      inibeta = rep(0,p)
      inibeta[order] = tempbeta
    }
    if(style=="classic"){
      temprss = RSS(X, y, inibeta)
      tempbeta = beta_EM(X=X,y=y,inibeta=inibeta,M=M,a=a)    
      deltabeta = tempbeta - inibeta
      testrss = c()
      for(k in 1:4){
        testbeta = k*deltabeta + inibeta
        testrss = c(testrss, RSS(X, y, testbeta))
      }
      k = which.min(testrss)
      tempbeta = k*deltabeta + inibeta
      order = S_M(tempbeta, M)
      inibeta = rep(0,p)
      inibeta[order] = tempbeta[order]
    }
    oldrss = rss
    rss=sum((y-X%*%as.matrix(inibeta))^2)
    if(abs(rss-oldrss)<10^(-4)||(oldrss<rss)){
      break
    }
    if(print==T){
      print(rss)
    }
  }
  time = proc.time() - time
  hatset = sort((1:p)[inibeta!=0])
  return(list(inibeta,tempbeta,rss=rss,time=time,hatset=hatset))    
}


combineradius = function(X, y, M,flag=200,print=F,inipoint,iniTem=20,annealing=F,radiusceil=10,cutrate=2){
  time = proc.time()
  p = dim(X)[[2]]
  n = dim(X)[[1]]
  tempvar = 1:p
  tempbeta = subset_EMtest(X, y, M, flag=1000,print=F,style="classic",ini="other",inipoint=inipoint)[[2]]
  temporder = order(abs(tempbeta),decreasing=T)
  temporder[1:(M)] = temporder[(M):1]  
  tmp = rep(0,p)
  tmp[tempvar] = tempbeta
  tempbeta = tmp
  tempvar = tempvar[temporder]
  
  castnum = floor((p - M - 1)/cutrate)
  tempvar = tempvar[1:(p-castnum)]
  
  kset = c(0)
  varlist = list(tempvar)
  betalist = list(tempbeta)
  
  hatbeta = gliner(X[,tempvar[1:M]],y)
  rss = RSS(X[,tempvar[1:M]],y,hatbeta)
  bestset = tempvar[1:M]
  rssset = c(rss)
  iter = 1
  stageset = c(1)
  BREAK = 0
  while(length(kset)>0){
    if(BREAK){
      break
    }
    fathervar = varlist[[1]]
    varlist = varlist[-1]
    fatherk = kset[1]
    kset = kset[-1]
    fatherbeta = betalist[[1]]
    betalist = betalist[-1]
    fatherrss = rssset[1]
    rssset = rssset[-1]
    stage = stageset[1]
    stageset = stageset[-1]
    fathern = length(fathervar)
    
    for(tempk in fatherk:(min(M,fathern-2))){
      tempvar = fathervar[-(tempk+1)]
      inibeta = fatherbeta[tempvar]
      Tem = iniTem * (0.4)^(stage-1)
      radius = stage + tempk
      ns = length(tempvar)
      if(tempk >= M){
        next
      }
      if(ns<M){
        BREAK = 1
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
        BREAK = 1
        break
      }
      if(radius<radiusceil){
        tempbeta = subset_EMtest(X[,tempvar], y, M, flag=iter,print=F,style="classic",ini="other",inipoint=inibeta)[[2]]
        castnum = floor((ns - M - 1)/cutrate)
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
        cat(flag,temprss,rss,stage,length(varlist),"\n")
      }
      if(ns==M){
        BREAK = 1
        break
      }
      if(tempk>(min(M,ns-2))){
        next
      }
      kset = c(kset, tempk)
      varlist = c(varlist, list(tempvar))
      betalist = c(betalist, list(tempbeta))
      rssset = c(rssset, temprss)
      stageset = c(stageset, stage+1)
    }
  }
  time = proc.time() - time
  return(list(hatset=bestset,rss=rss,time=time))
}










