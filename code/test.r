 

set2 = c(sample((1:100)[-trueset],10), trueset)
set1 = c(sample((1:100)[-trueset],36), trueset[-1])
beta1 = gliner(X[,set1],y)
beta2 = gliner(X[,set2],y)
RSS(X[,set1],y,beta1)
RSS(X[,set2],y,beta2)

n = 5
p = 5
truelength = 2
M = 3
X = generate(n,p,0.5)
trueset = sample(1:p,truelength)
print(sort(trueset))
truebeta = rep(0,p)
truebeta[trueset] = rep(3,truelength)
epsilon = rnorm(n,0,1.5)
y = X %*% as.matrix(truebeta) + epsilon 

combine(X,y,M)






