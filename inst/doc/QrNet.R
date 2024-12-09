## ----eval=FALSE---------------------------------------------------------------
#  svda <- irlba::irlba(A, nu = r, nv = r)
#  uu <- svda$u[,1]
#  vv <- svda$v[,1]
#  d <- svda$d[1]

## ----eval=FALSE---------------------------------------------------------------
#  ret_two_stage <- lm.fit(cbind(X, uu[1:n_obs], vv[1:n_obs]), y)
#  ret_two_stage$beta <- ret_two_stage$coefficients

## ----include=FALSE------------------------------------------------------------
vec_norm <- function(x) {sqrt(sum(x^2))}   

## -----------------------------------------------------------------------------
library(SA24204148)
n <- 100
p <- 3
sigmaa <- 1
sigmay <- 1e-5
A <- matrix(rnorm(n^2, sd = sigmaa), nrow = n)
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
y <- rnorm(n, sd = sigmay)
ret <- two_stage(A, X, y)
ret$coefficients
vec_norm(ret$residuals)

## ----eval=FALSE---------------------------------------------------------------
#  int update_beta(const arma::colvec& y,
#                  const arma::mat& X,
#                  const arma::colvec& u,
#                  const arma::colvec& v,
#                  arma::colvec& beta,
#                  const arma::colvec& theta,
#                  double sigma,
#                  const arma::colvec& z){
#    int n_obs = X.n_rows;
#    arma::mat WW = arma::join_rows(X, u.head(n_obs));
#    arma::mat W = arma::join_rows(WW, v.head(n_obs));
#    beta = arma::solve(W, y - z + theta / sigma);
#    return 1;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  int update_z(const arma::colvec& y,
#               const arma::mat& X,
#               const arma::colvec& u,
#               const arma::colvec& v,
#               const arma::colvec& beta,
#               const arma::colvec& theta,
#               double sigma,
#               double tau,
#               arma::colvec& z){
#    int n_obs = X.n_rows;
#    arma::mat W = join_rows(X, u.head(n_obs), v.head(n_obs));
#    arma::colvec xi = y -  W * beta + theta/sigma;
#    for(int i = 0; i < n_obs; i++){
#      z(i) = prox(tau, xi(i), n_obs * sigma);
#    }
#    return 1;
#  };

## ----eval=FALSE---------------------------------------------------------------
#  double update_d(const arma::mat& A,
#                  const arma::colvec& u,
#                  const arma::colvec& v)
#  {
#    double d = as_scalar(u.t() * A * v);
#    return d;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  int update_u(const arma::mat& A,
#               const arma::colvec& y,
#               const arma::mat& X,
#               arma::colvec& u,
#               const arma::colvec& v,
#               const arma::colvec& beta,
#               const arma::colvec& theta,
#               double d,
#               double sigma,
#               double tau,
#               double l,
#               const arma::colvec& z){
#    int n_obs = X.n_rows, k = X.n_cols, n = A.n_rows, n_c = n - n_obs;
#    arma::colvec rmat = beta(k)*(sigma*(y - X*beta.head(k) - beta(k+1)*v.head(n_obs) - z) + theta) + 2*l*d*A.head_rows(n_obs)*v/n;
#    u.head(n_obs) = rmat / (pow(beta(k),2)*sigma + 2*l*pow(d,2));
#  
#    if(n_obs < n) {
#      u.tail(n_c) = A.tail_rows(n_c)*v/d/n;
#    }
#  
#    // normalize u
#    u = arma::normalise(u) * sqrt(n);
#  
#    return 1;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  int update_v(const arma::mat& A,
#               const arma::colvec& y,
#               const arma::mat& X,
#               const arma::colvec& u,
#               arma::colvec& v,
#               const arma::colvec& beta,
#               const arma::colvec& theta,
#               double d,
#               double sigma,
#               double tau,
#               double l,
#               const arma::colvec& z){
#    int n_obs = X.n_rows, k = X.n_cols, n = A.n_rows, n_c = n - n_obs;
#    arma::colvec rmat = beta(k+1)*(sigma*(y - X*beta.head(k) - beta(k)*u.head(n_obs) - z) + theta) + 2*l*d*A.head_cols(n_obs).t()*u/n;
#    v.head(n_obs) = rmat / (pow(beta(k+1),2)*sigma + 2*l*pow(d,2));
#  
#    if(n_obs < n) {
#      v.tail(n_c) = A.tail_cols(n_c).t()*u/d/n;
#    }
#  
#    // normalize u
#    v = arma::normalise(v) * sqrt(n);
#  
#    return 1;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  int update_theta(arma::colvec& theta,
#                   double sigma,
#                   const arma::mat& X,
#                   const arma::colvec& u,
#                   const arma::colvec& v,
#                   const arma::colvec& y,
#                   const arma::colvec& z,
#                   const arma::colvec& beta){
#    int n_obs = X.n_rows;
#    arma::mat WW = arma::join_rows(X, u.head(n_obs));
#    arma::mat W = arma::join_rows(WW, v.head(n_obs));
#    theta -= sigma*(W*beta+z-y);
#    return 1;
#  }

## -----------------------------------------------------------------------------
library(irlba)
n <- 1000
p <- 3
sigmaa <- 1
sigmay <- 1e-5
A <- matrix(rnorm(n^2, sd = sigmaa), nrow = n)
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
y <- rnorm(n, sd = sigmay)
ret_q <- QrNet(A, X, y)
ret_q$beta

## -----------------------------------------------------------------------------
set.seed(0)
n=1e3
d=1
p=3

sigma_a_set=c(2^(-4),2^(-2))
u=matrix(rnorm(n,0,1),ncol = 1)
v=0.5*u+matrix(rnorm(n,0,1),ncol = 1)

u = u/sqrt(sum(u^2))*sqrt(n)
v = v/sqrt(sum(v^2))*sqrt(n)

sigma_a=sigma_a_set[2] #  choose sigma_a
A=d*u%*%t(v)+matrix(rnorm(n^2,0,sigma_a^2),n,n)

beta_x=matrix(c(1,3,5),ncol=1)
beta_v=1
beta_u_set=c(2^0,2^2,2^4)
sigma_y_set=c(2^(-4),2^(-2),2^0)
beta_u=beta_u_set[1] #  choose beta_u
sigma_y=sigma_y_set[3] #  choose beta_y

X=matrix(rnorm(n*p,0,1),nrow = n)
y=X%*%beta_x+u*beta_u+v*beta_v+matrix(rnorm(n,0,sigma_y^2),ncol = 1)

lambda0=n*sigma_y^2/sigma_a^2

# 1. 查看生成的目标变量 y 和预测变量 X 的关系

# 可视化 X 和 y 的关系（散点图）
plot(y, type = "p", col = "blue", pch = 16, main = "Target variable y vs. Index", xlab = "Index", ylab = "y")

# 可视化 X 中的前三个变量
matplot(X[, 1:3], type = "l", col = c("red", "green", "blue"), lty = 1, 
        main = "First Three Predictors (X)", xlab = "Index", ylab = "Value")
legend("topright", legend = c("X1", "X2", "X3"), col = c("red", "green", "blue"), lty = 1)

# 2. 可视化网络数据 A 的结构（前几行）
image(A[1:20, 1:20], main = "Network Structure A (Subset)", xlab = "Index", ylab = "Index")

# 3. 可视化 u 和 v 的关系

# u 和 v 的关系（散点图）
plot(u, v, col = "purple", pch = 16, main = "Relationship between u and v", 
     xlab = "u", ylab = "v")

# u 和 v 的时间序列图
plot(u, type = "l", col = "red", main = "Time Series of u and v", xlab = "Index", ylab = "Value")
lines(v, col = "blue")
legend("topright", legend = c("u", "v"), col = c("red", "blue"), lty = 1)

## -----------------------------------------------------------------------------
library(knitr)  # 用于生成表格

# 运行 two_stage 方法
ret <- two_stage(A, X, y)
beta_two_stage <- ret$beta

# 运行 QrNet 方法
tau_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
beta_qrnet <- lapply(tau_values, function(tau) {
  ret_q <- QrNet(A, X, y, tau = tau)
  ret_q$beta
})

# 将结果整理为数据框
result_table <- data.frame(
  Beta = paste0("Beta_", 1:length(beta_two_stage)),
  Two_Stage = beta_two_stage
)

# 添加不同分位数的 QrNet 结果到数据框
for (i in seq_along(tau_values)) {
  result_table[[paste0("Tau = ", tau_values[i])]] <- beta_qrnet[[i]]
}

# 打印表格
kable(result_table, format = "markdown", caption = "Comparison of Two-Stage and QrNet Methods")

