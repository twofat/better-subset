rm(list=ls())
setwd("E://Variable Selective Methods//better subset")
source("E://Variable Selective Methods//better subset//functions.r")

n=100
p=300
rho = 0.7
truelength = 10
M = 20
sd = 1
random = randomdesign(n,p,rho,truelength,sd)
X = random$X
y = random$y
trueset = random$trueset

time = proc.time()
betahat = subset_EM(X, y, M, flag=10000,print=T,style="classic",ini="zero")[[1]]
hatset1 = sort((1:p)[betahat!=0])
hatset1
sum(trueset %in% hatset1 / truelength)
print(proc.time()-time)

time = proc.time()
betahat = subset_EM(X, y, M, flag=10000,print=T,style="classic",ini="")[[1]]
hatset3 = sort((1:p)[betahat!=0])
hatset3
sum(trueset %in% hatset1 / truelength)
print(proc.time()-time)

source("E://Variable Selective Methods//better subset//functions.r")
time = proc.time()
hatset2 = sort(combineradius(X, y, M,flag=200,print=T,inipoint=betahat,annealing=T,iniTem=5,radiusceil=10)[[1]])
print(proc.time()-time)

time = proc.time()
hatset2 = sort(combinetest(X, y, M,flag=100,print=T,inipoint=betahat,annealing=F,Tem=20)[[1]])
hatset2 = sort(combineradius(X, y, M,flag=100,print=T,inipoint=betahat,annealing=F,radiusceil=1000,iniTem=20)[[1]])
print(proc.time()-time)

combinesearch(X, y, M,flag=500,print=T,inipoint=betahat,annealing=F,Tem=20)[[1]]

set = c(trueset, sample((1:p)[-trueset],M-length(trueset)))
hatset =set
hatbeta = gliner(X[,hatset],y)
RSS(X[,hatset],y,hatbeta)




rm(list=ls())

source("E://Variable Selective Methods//better subset//iterative_test.r")
RSS = matrix(apply(rss, 2, mean),,3,byrow=T)
RATE = matrix(apply(100*correctrate, 2, mean),,3,byrow=T)
rownames(RSS) = methods1
rownames(RATE) = methods1
colnames(RSS) = c("ini", "EM-Impr", "H-Radius")
colnames(RATE) = c("ini", "EM-Impr", "H-Radius")
RATE = round(RATE,1) 
RSS
RATE


save(RSS,RATE,rss,correctrate,file="simulate200500sd15.RData")

load(file="simulate50100.RData")


