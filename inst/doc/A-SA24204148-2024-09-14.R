## ----include=TRUE,echo=FALSE--------------------------------------------------
library("latex2exp")
library("ggplot2")

## -----------------------------------------------------------------------------
rRayleigh <- function(n = 1000, sigma, seed = 0){
  f <- function(x, sigma){
    return(x/sigma^2*exp(-x^2/2/sigma^2))
  }
  F_inverse <- function(x, sigma){
    return(sqrt(-2*sigma^2*log(1-x))) 
  }
  set.seed(seed)
  u <- runif(n)
  x <- sapply(u,function(t) F_inverse(t,sigma))
  return(x)
}

n = 5000
x_1 = rRayleigh(n = n, sigma = 1)
p = sqrt(-2*log(1-seq(0.1,0.9,0.2)))
knitr::kable(round(rbind(quantile(x_1,probs = seq(0.1,0.9,0.2)),p),3))


## -----------------------------------------------------------------------------
op <- par(mfrow = c(2, 3),pty = "s")       
sigma_set = c(0.5,1,2,5,10,20)
for(sigma in sigma_set){
  x = rRayleigh(n = 5000, sigma = sigma)
  hist(x,breaks = 30,prob = T, main = bquote(sigma == ~ .(sigma)), xlab = "x",col = "white")
  y = seq(0,max(x),.01)
  lines(y,sapply(y,function(t) t/sigma^2*exp(-t^2/2/sigma^2)),col="red")
}
par(op)

## -----------------------------------------------------------------------------
normal_mixture <- function(n,p1,seed=0){
  set.seed(0)
  z = sample(0:1, size = n, replace = T, prob = c(p1,1-p1))
  x = rnorm(n)+z*3
  return(x)
}

n = 1000
x = normal_mixture(n, 0.75)
plot(density(x),ylim=c(0,0.4),main = TeX('$p_1=0.75$'),xlab = "",col="blue")
hist(x,breaks = 50,add=T,probability = T,col = "white")
curve(1/sqrt(2*pi)*exp(-x^2/2),add = T,col="gray")
curve(1/sqrt(2*pi)*exp(-(x-3)^2/2),add = T,col="gray")

## -----------------------------------------------------------------------------
plot(density(normal_mixture(1000, 0.5)),ylim=c(0,0.4),main = "",xlab = "x")
for(p in c(seq(0.01,0.1,0.02), seq(0.1,0.4,0.1))){
  lines(density(normal_mixture(1000, p)))
}

## -----------------------------------------------------------------------------
cpp <- function(n,t,lambda,alpha,beta,seed=0){
  set.seed(seed)
  x = sapply(rpois(n,lambda*t),function(t) sum(rgamma(t,alpha,beta)))
  return(x)
}

## -----------------------------------------------------------------------------
op <- par(mfrow = c(1, 2),pty = "s")
t = 10
n = 1000
lambda = 1
alpha = c(1,2,3,4)
beta = seq(0.5,10,0.5)
for (i in seq(length(alpha))) {
  x = sapply(beta,function(z) cpp(n,t,lambda,alpha[i],z))
  if(i==1){plot(x=beta,y=apply(x,2,mean),col=i,xlab=bquote(beta),ylab="mean",pch=2)}
  else {points(beta,y=apply(x,2,mean),col=i,pch=2)}
  curve(i/x*lambda*t,add=T,col=i)
}
legend("topright",legend = sapply(seq(4),function(t) bquote(alpha==~.(t))),col = seq(4),lty=c(1,1,1,1))

for (i in seq(length(alpha))) {
  x = sapply(beta,function(z) cpp(n,t,lambda,alpha[i],z))
  if(i==1){plot(x=beta,y=apply(x,2,var),col=i,xlab=bquote(beta),ylab="var",pch=5,cex=0.5)}
  else {points(beta,y=apply(x,2,var),col=i,pch=5,cex=0.5)}
  curve(i*(1+i)/x^2*lambda*t,add=T,col=i)
}
legend("topright",legend = sapply(seq(4),function(t) bquote(alpha==~.(t))),col = seq(4),lty=c(1,1,1,1))


## -----------------------------------------------------------------------------
alpha = 1
beta = 1/2
lambda = c(1,3,8,20)
n = c(10,100,1000,10000)
mean_hat = NULL
var_hat = NULL
for(j in seq(length(n))){
  x = sapply(lambda,function(z) cpp(n[j],t,z,alpha,beta))
  mean_hat = rbind(mean_hat,apply(x,2,mean))
  var_hat = rbind(var_hat,apply(x,2,var))
}

mean_hat = rbind(mean_hat,lambda*t*alpha/beta)
colnames(mean_hat) <- lambda
rownames(mean_hat) <- c(n,"True")
mean_hat
var_hat = rbind(var_hat,lambda*t*alpha*(1+alpha)/beta^2)
colnames(var_hat) <- lambda
rownames(var_hat) <- c(n,"True")
var_hat

