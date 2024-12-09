## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
f0 <- function(x) return(1)
f1 <- function(x) return(exp(-x))
f2 <- function(x) 1/pi/(1+x^2)
f3 <- function(x) return(1/exp(x)/(1-1/exp(1)))
f4 <- function(x) 4/pi/(1+x^2)
g<- function(x){
  ##  exp(-x - log(1+x^2)) * (x > 0) * (x < 1) 
  if(x<1 && x>0){
    return(1/(1+x^2)/exp(x))
  }
  else return(0)
}
curve(f1,xlim = c(0,1),ylim=c(0,2),lty=3,col="green",ylab="")
curve(f2,add=T,lty=4,col="skyblue")
curve(f3,add=T,lty=5,col="orange")
curve(f4,add=T,lty=6,col="red")
x = seq(0.01,1-0.01,0.01)
points(x,sapply(x, g),type = "l")
points(x,sapply(x, f0),type="l",lty=2,col="gray")
legend("topright",legend = c("g",as.character(0:4)),lty=c(1:6),col = c(1,"gray","green","skyblue","orange","red"),ncol=2)

## -----------------------------------------------------------------------------
set.seed(0)
m = 10000

g<- function(x){
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1) 
}

theta.hat <- se <- numeric(5) 
x <- runif(m) #using f0
fg <- g(x) 
theta.hat[1] <- mean(fg) 
se[1] <- sd(fg)
x <- rexp(m, 1) #using f1
fg <- g(x) / exp(-x) 
theta.hat[2] <- mean(fg) 
se[2] <- sd(fg)
x <- rcauchy(m) #using f2 
i <- c(which(x > 1), which(x < 0)) 
x[i] <- 2 #to catch overflow errors in g(x) 
fg <- g(x) / dcauchy(x) 
theta.hat[3] <- mean(fg) 
se[3] <- sd(fg)
u <- runif(m) #f3, inverse transform method 
x <- - log(1 - u * (1 - exp(-1))) 
fg <- g(x) / (exp(-x) / (1 - exp(-1))) 
theta.hat[4] <- mean(fg) 
se[4] <- sd(fg)
u <- runif(m) #f5, inverse transform method 
x <- tan(pi * u / 4) 
fg <- g(x) / (4 / ((1 + x^2) * pi)) 
theta.hat[5] <- mean(fg) 
se[5] <- sd(fg)
knitr::kable(rbind(theta.hat, se))

## -----------------------------------------------------------------------------
data0 = matrix(c(24,35,21,30,1355,603,192,224),byrow = F,ncol = 2)
data=cbind(data0,data0[,1]/(data0[,1]+data0[,2]))

s2 = c(0,1,2,3)
s1 = s2*2
s3 = s2+1

x = rep(rep(s1,2),as.vector(data0))
y = rep(c(rep(1,4),rep(0,4)),as.vector(data0))
fit.lm_a = lm(y~x)
fit.log_a<-glm(data0 ~ s1, family=binomial)

x = rep(rep(s2,2),as.vector(data0))
y = rep(c(rep(1,4),rep(0,4)),as.vector(data0))
fit.lm_b = lm(y~x)
fit.log_b<-glm(data0 ~ s2, family=binomial)

x = rep(rep(s3,2),as.vector(data0))
y = rep(c(rep(1,4),rep(0,4)),as.vector(data0))
fit.lm_c = lm(y~x)
fit.log_c<-glm(data0 ~ s3, family=binomial)
        
plot(s1,data[,3],ylim=c(0,0.25),col=2,xlab='Snoring',ylab='prop',pch=16)
points(s2,data[,3],col=3,pch=16)
points(s3,data[,3],col=4,pch=16)
xs<-seq(0,6,0.01)
pr.lm_a<-fit.lm_a$coef[1]+xs*fit.lm_a$coef[2]
lines(xs,pr.lm_a,col=2,lwd=2)
pr.log_a<-fit.log_a$coef[1] + fit.log_a$coef[2]*xs
lines(xs,exp(pr.log_a)/(1+exp(pr.log_a)),col=2,lwd=2,lty=2)
pr.lm_b<-fit.lm_b$coef[1]+xs*fit.lm_b$coef[2]
lines(xs,pr.lm_b,col=3,lwd=2)
pr.log_b<-fit.log_b$coef[1] + fit.log_b$coef[2]*xs
lines(xs,exp(pr.log_b)/(1+exp(pr.log_b)),col=3,lwd=2,lty=2)
pr.lm_c<-fit.lm_c$coef[1]+xs*fit.lm_c$coef[2]
lines(xs,pr.lm_c,col=4,lwd=2)
pr.log_c<-fit.log_c$coef[1] + fit.log_c$coef[2]*xs
lines(xs,exp(pr.log_c)/(1+exp(pr.log_c)),col=4,lwd=2,lty=2)
legend("topleft",legend = c("(a)linear","(b)linear","(c)linear","(a)logistic","(b)logistic","(c)logistic"),ncol=2,lty=c(1,1,1,2,2,2),col=c(2,3,4),cex = 0.7)

## -----------------------------------------------------------------------------
beta = matrix(c(fit.lm_a$coefficients[2],fit.lm_b$coefficients[2],fit.lm_c$coefficients[2],fit.log_a$coefficients[2],fit.log_b$coefficients[2],fit.log_c$coefficients[2]),nrow=2,byrow = T)
colnames(beta)=c("a","b","c")
rownames(beta)=c("linear","logistic")
knitr::kable(beta)

## -----------------------------------------------------------------------------
set.seed(0)
n = 10000
Y = runif(n)*2*pi
X1 = cos(Y)
hist(X1,breaks = 50)

## -----------------------------------------------------------------------------
f <- function(x) 1/pi/sqrt(1-x^2)
hist(X1,breaks = 50,freq = F)
curve(f,from = -1,to = 1,add=T,col="red")

