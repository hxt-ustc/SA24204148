## -----------------------------------------------------------------------------
# 计算样本偏度
skewness <- function(sample){
  s = sqrt(var(sample))
  x_bar = mean(sample)
  b1 = mean(((sample-x_bar)/s)^3)
  return(b1)
}

# 设置参数
N_sim <- 1000  # 模拟次数
n <- 100  # 样本大小

# 初始化用于存储偏度的向量
skewness_values <- numeric(N_sim)

# 进行蒙特卡洛模拟
set.seed(123)  # 设置随机种子以便结果可重复
for (i in 1:N_sim) {
  # 从标准正态分布中生成样本
  sample <- rnorm(n)
  
  # 计算偏度（b_1 的平方根）
  skewness_values[i] <- skewness(sample) #sqrt(skewness(sample))
}

# 估计所需的分位数
quantiles <- quantile(skewness_values, probs = c(0.025, 0.05, 0.95, 0.975))

# 计算分位数估计的标准误差
# 使用近似公式 Var(hat{x}_q) = (q * (1 - q)) / (n * f(x_q)^2)
# f(x_q) 在此用正态分布的密度近似
se_estimates <- numeric(length(quantiles))

for (i in seq_along(quantiles)) {
  q <- c(0.025, 0.05, 0.95, 0.975)[i]
  f_xq <- dnorm(quantiles[i], mean = 0, sd = 1)
  se_estimates[i] <- sqrt(q * (1 - q) / (N_sim * f_xq^2))
}

# 打印分位数估计及其标准误差
for (i in seq_along(quantiles)) {
  cat(sprintf("Quantile %f: Estimate = %f, SE = %f\n", 
              c(0.025, 0.05, 0.95, 0.975)[i], quantiles[i], se_estimates[i]))
}

# 比较估计的分位数和大样本近似的分位数
approx_quantiles <- qnorm(c(0.025, 0.05, 0.95, 0.975), mean = 0, sd = sqrt(6 / n))
cat("\nLarge Sample Approximation Quantiles:\n")
print(approx_quantiles)

result = cbind(quantiles,approx_quantiles,se_estimates)
colnames(result) <- c("estimate", "approx", "se")
knitr::kable(round(result,3))

## -----------------------------------------------------------------------------
# 设置参数
N_sim <- 10000  # 模拟次数
n <- 50  # 样本大小
alpha <- 0.05  # 显著性水平

# 初始化用于存储 p 值的矩阵
p_values <- matrix(0, nrow = N_sim, ncol = 3)
colnames(p_values) <- c("Pearson", "Spearman", "Kendall")

set.seed(123)  # 设置随机种子

# 1. 双变量正态分布下的功效比较
for (i in 1:N_sim) {
  # 生成双变量正态分布样本
  X <- rnorm(n)
  Y <- 0.5 * X + sqrt(1 - 0.5^2) * rnorm(n)
  
  # 使用 cor.test() 进行三种相关性检验
  p_values[i, "Pearson"] <- cor.test(X, Y, method = "pearson")$p.value
  p_values[i, "Spearman"] <- cor.test(X, Y, method = "spearman")$p.value
  p_values[i, "Kendall"] <- cor.test(X, Y, method = "kendall")$p.value
}

# 计算功效（拒绝原假设的比例）
power_pearson <- mean(p_values[, "Pearson"] < alpha)
power_spearman <- mean(p_values[, "Spearman"] < alpha)
power_kendall <- mean(p_values[, "Kendall"] < alpha)

cat("双变量正态分布下的经验功效:\n")
cat(sprintf("Pearson: %f, Spearman: %f, Kendall: %f\n", power_pearson, power_spearman, power_kendall))

## -----------------------------------------------------------------------------
# 2. 非正态分布下非参数检验的高功效示例
# 生成一个非正态分布的双变量依赖数据集
set.seed(123)
p_values_non_normal <- matrix(0, nrow = N_sim, ncol = 3)
colnames(p_values_non_normal) <- c("Pearson", "Spearman", "Kendall")

for (i in 1:N_sim) {
  # 生成具有依赖关系的非正态分布样本，例如 Y = X^2 + 噪声
  X <- runif(n, -1, 1)
  Y <- X^2 + rnorm(n, mean = 0, sd = 0.1)
  
  # 使用 cor.test() 进行三种相关性检验
  p_values_non_normal[i, "Pearson"] <- cor.test(X, Y, method = "pearson")$p.value
  p_values_non_normal[i, "Spearman"] <- cor.test(X, Y, method = "spearman")$p.value
  p_values_non_normal[i, "Kendall"] <- cor.test(X, Y, method = "kendall")$p.value
}

# 计算非正态分布下的功效
power_pearson_non_normal <- mean(p_values_non_normal[, "Pearson"] < alpha)
power_spearman_non_normal <- mean(p_values_non_normal[, "Spearman"] < alpha)
power_kendall_non_normal <- mean(p_values_non_normal[, "Kendall"] < alpha)

cat("\n非正态分布下的经验功效:\n")
cat(sprintf("Pearson: %f, Spearman: %f, Kendall: %f\n", power_pearson_non_normal, power_spearman_non_normal, power_kendall_non_normal))

