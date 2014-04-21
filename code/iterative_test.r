setwd("E://Variable Selective Methods//better subset")
source("E://Variable Selective Methods//better subset//functions.r")

n = 200
p = 500
rho = 0.5
truelength = 10
sd = 1.5
M = 20

rss = c()
time = c()
correctrate = c()
methods1 = c("SIS","SF","Lar","LASSO","Elastic","GLM","SCAD")
inimethods = c("zero", "stepforward", "lar", "lasso", "elastic","glm","scad")
methods = c()
for(method in methods1){
  methods = c(methods, method)
  methods = c(methods, paste("Im", method, sep=""))
  #  methods = c(methods, paste("BBA", method, sep=""))
  #  methods = c(methods, paste("RBBA", method, sep=""))
  methods = c(methods, paste("HRBBA", method, sep=""))
}

for(iterativenum in 1:100){
  cat("iterativenum:",iterativenum,"\n")
  ran = randomdesign(n=n, p=p, rho=rho, truelength=truelength, sd=sd)
  X = ran$X
  y = ran$y
  trueset = ran$trueset
  thisresult = list()
  for(ini in inimethods){
    ori = subset_EM(X, y, M, flag=0,print=F,style="classic",ini=ini)
    Im = subset_EM(X, y, M, flag=10000,print=F,style="classic",ini=ini)
#    BBA = combineradius(X, y, M,flag=500,print=F,inipoint=Im[[1]],annealing=F,iniTem=5,radiusceil=1000)
#    RBBA = combineradius(X, y, M,flag=500,print=F,inipoint=Im[[1]],annealing=F,iniTem=5,radiusceil=10)
    HRBBA = combineradius(X, y, M,flag=500,print=F,inipoint=Im[[1]],annealing=T,iniTem=5,radiusceil=10) 
    thisresult = c(thisresult, list(ori, Im, HRBBA))
  }
  
  tempcorrectrate = c()
  temptime = c()
  temprss = c()
  tempcorrectrate = c()
  for(k in thisresult){
    temptime = c(temptime, k$time[[1]])
    temprss = c(temprss, k$rss)
    tempcr = sum(trueset %in% (k$hatset)) / truelength
    tempcorrectrate = c(tempcorrectrate, tempcr)
  }
  tempprint = matrix(tempcorrectrate,,3,byrow=T)
  rownames(tempprint) = methods1
  print(tempprint)
  rss = rbind(rss, temprss)
  time = rbind(time, temptime)
  correctrate = rbind(correctrate, tempcorrectrate)
}


colnames(rss) = methods
colnames(time) = methods
colnames(correctrate) = methods









