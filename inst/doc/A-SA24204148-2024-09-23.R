## -----------------------------------------------------------------------------
beta33_cdf <- function(x, n =1000, seed=0){
  set.seed(seed)
  u = runif(n, max = x)
  return(x*mean(30*u^2*(1-u)^2))
}

## -----------------------------------------------------------------------------
t = seq(0.1,0.9,0.1)
cdf.hat = sapply(t,beta33_cdf)
true.cdf = pbeta(t,3,3)
result = rbind(cdf.hat,true.cdf)
colnames(result) = t
result = rbind(result,cdf.hat-true.cdf)
result

## -----------------------------------------------------------------------------
rRayleigh <- function(sigma,n=1200,seed=0){
  set.seed(seed)
  u = runif(n)
  return(list(x1 = sqrt(-2*sigma^2*log(1-u)), x2 = sqrt(-2*sigma^2*log(u))))
}

## -----------------------------------------------------------------------------
sigma = c(1,2,5,10,20)
f1 <- function(sigma,n=1200,seed=0){
  x1 = rRayleigh(sigma,n,seed)
  x1 = (x1$x1 + x1$x2)/2
  x2 = (rRayleigh(sigma,n,seed)$x1 + rRayleigh(sigma,n,2*seed)$x1)/2
  return(var(x1)/var(x2))
}
round(1-sapply(sigma,f1),3)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))  # 设置图形区域为1行2列

sigma = 1

X_plus_independent = (rRayleigh(sigma, seed = 1)$x1 + rRayleigh(sigma,seed = 2)$x1)/2
X_plus_antithetic = rRayleigh(sigma)
X_plus_antithetic = (X_plus_antithetic$x1 + X_plus_antithetic$x2)/2
    
hist(X_plus_independent, breaks = 30, col = 'deepskyblue', 
         main = '独立变量生成的Rayleigh样本', 
         xlab = '样本值', 
         ylab = '频率', 
         xlim = c(0, max(X_plus_independent, X_plus_antithetic)),
         probability = TRUE)
    
hist(X_plus_antithetic, breaks = 30, col = 'deeppink', 
         main = '对立变量生成的Rayleigh样本', 
         xlab = '样本值', 
         ylab = '频率', 
         xlim = c(0, max(X_plus_independent, X_plus_antithetic)),
         probability = TRUE)
    
legend("topright", legend = c("独立变量", "对立变量"), 
           fill = c("deepskyblue", "deeppink"))

## -----------------------------------------------------------------------------
g <- function(x){
  return(x^2/sqrt(2*pi)*exp(-x^2/2))
}

f1 <- function(x){
  return(x*sqrt(exp(1))*exp(-(x)^2/2))
}

f2 <- function(x){
  return(dnorm(x)/(1-pnorm(1)))
}

curve(g,xlim = c(1,3),ylim=c(0,2),ylab = "")
curve(f1,add=T,col="red")
curve(f2,add=T,col="green")
legend("topright",legend = c("g","f1","f2"),col = c(1,"red","green"),lty=1)

## -----------------------------------------------------------------------------
F1 <- function(x){
  return(sqrt(1-2*log(1-x)))
}

F2 <- function(x){
  return(qnorm(pnorm(1)+x*(1-pnorm(1))))
}

n = 1000
set.seed(0)
u = runif(n)
x1 = F1(u)
x2 = F2(u)
y1 = g(x1)/f1(x1)
y2 = g(x2)/f2(x2)
mean(y1)
mean(y2)
var(y1)
var(y2)

## -----------------------------------------------------------------------------
fast_sort <- function(x){
  if(length(x)<2) return(x)
  x0 <- sample(x,1)
  return(c(fast_sort(x[x<x0]),x0,fast_sort(x[x>x0])))
}

## -----------------------------------------------------------------------------
n = c(1,2,4,6,8)*1e4
simulation = 150
times <- function(n,simulation){
  t = 0
  for(i in seq(simulation))  {
    x = sample(seq(n),n)
    start_time <- Sys.time()
    fast_sort(x)
    t = t + Sys.time()-start_time
    }
  return(t/simulation)
}

t = sapply(n,function(z) times(z,simulation))

## -----------------------------------------------------------------------------
y = n*log(n)
my.fit <- lm(t~y)
plot(n*log(n),t,ylab=expression(an))
curve(my.fit$coefficients[1]+my.fit$coefficients[2]*x,add = T)

