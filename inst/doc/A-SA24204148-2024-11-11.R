## ----include=FALSE------------------------------------------------------------
library(lpSolve)
library(microbenchmark)

## -----------------------------------------------------------------------------
# 设置目标函数系数（最小化 4x + 2y + 9z）
objective <- c(4, 2, 9)

# 设置约束条件系数
constraints <- matrix(c(2, 1, 1,   # 2x + y + z <= 2
                        1, -1, 3), # x - y + 3z <= 3
                      nrow = 2, byrow = TRUE)

# 设置约束条件的右侧值
rhs <- c(2, 3)

# 设置约束类型
direction <- c("<=", "<=")

# 使用 lpSolve 求解
solution <- lp("min", objective, constraints, direction, rhs, all.int = FALSE, all.bin = FALSE)

# 查看结果
solution$objval      # 最小值
solution$solution    # x, y, z 的取值

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

models_for <- list()


for (i in seq_along(formulas)) {
  models_for[[i]] <- lm(formulas[[i]], data = mtcars)
}

lapply(models_for, summary)

## -----------------------------------------------------------------------------
models_lapply <- lapply(formulas, function(f) lm(f, data = mtcars))

lapply(models_lapply, summary)

## -----------------------------------------------------------------------------
index <- cbind(lapply(models_for, AIC),lapply(models_for, BIC))
rownames(index) <- as.character(formulas)
index

## -----------------------------------------------------------------------------
plot(lm(formulas[[4]],mtcars))

## -----------------------------------------------------------------------------
set.seed(0)
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

models_for <- list()

for (i in seq_along(bootstraps)) {
  models_for[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}

models_for

## -----------------------------------------------------------------------------
set.seed(0)
generate_bootstrap_sample <- function(data) {
  rows <- sample(1:nrow(data), rep = TRUE)
  data[rows, ]
}

# 创建一个空列表用于存储 bootstrap 样本
bootstraps <- list()

# 使用 for 循环生成 10 个 bootstrap 样本
for (i in 1:10) {
  bootstraps[[i]] <- generate_bootstrap_sample(mtcars)
}


fit_model <- function(data) {
  lm(mpg ~ disp, data = data)
}

models_lapply_2 <- lapply(bootstraps, fit_model)

models_lapply_2

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

# 使用 lapply 提取每个模型的 R^2 值
rsq_values_lapply <- lapply(models_lapply, rsq)

# 将结果转换为数值向量（如果需要）
rsq_values_lapply <- unlist(rsq_values_lapply)

# 查看 R^2 值
cat("R2 for the four models fitted in Exercise 3:", rsq_values_lapply, "\n")

# 使用 lapply 提取每个模型的 R^2 值
rsq_values_lapply_2 <- lapply(models_lapply_2, rsq)

# 将结果转换为数值向量（如果需要）
rsq_values_lapply_2 <- unlist(rsq_values_lapply_2)

# 查看 R^2 值
cat("R2 for the ten models fitted in Exercise 4:", rsq_values_lapply_2, "\n")

## -----------------------------------------------------------------------------
set.seed(123)
# Simulate the trials
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

# Use sapply() with an anonymous function to extract the p-value from each trial
p_values <- sapply(trials, function(test) test$p.value)

# View the extracted p-values
hist(p_values, breaks = 30, freq = F)

## -----------------------------------------------------------------------------
par_lapply <- function(FUN, ..., FUN.VALUE, USE.NAMES = TRUE) {
  # Use Map() to apply FUN in parallel over all inputs
  results <- Map(FUN, ...)
  
  # Use vapply to simplify the results into a vector or matrix
  vapply(results, identity, FUN.VALUE, USE.NAMES = USE.NAMES)
}

## -----------------------------------------------------------------------------
# Define a function to sum pairs of numbers
sum_fun <- function(x, y) x + y

# Call par_lapply with two parallel lists
result <- par_lapply(sum_fun, list(1, 2, 3), list(4, 5, 6), FUN.VALUE = numeric(1))
print(result)  # Should output: [1] 5 7 9

## -----------------------------------------------------------------------------
fast_chisq_test <- function(x, y) {
  # Check that inputs are numeric vectors of the same length
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y)) {
    stop("Both inputs must be numeric vectors of the same length.")
  }
  
  # Create a contingency table
  contingency_table <- table(x, y)
  
  # Calculate row and column sums
  row_totals <- rowSums(contingency_table)
  col_totals <- colSums(contingency_table)
  total <- sum(contingency_table)
  
  # Compute expected values
  expected <- outer(row_totals, col_totals) / total
  
  # Compute chi-square statistic
  observed <- as.vector(contingency_table)
  expected <- as.vector(expected)
  chisq_stat <- sum((observed - expected)^2 / expected)
  
  return(chisq_stat)
}

## -----------------------------------------------------------------------------
# Load the microbenchmark package
library(microbenchmark)

# Define the fast version of chi-square test (as described earlier)
fast_chisq_test <- function(x, y) {
  if (length(x) != length(y)) {
    stop("The two vectors must be of the same length.")
  }
  
  # Calculate observed and expected values
  observed <- table(x, y)
  row_totals <- rowSums(observed)
  col_totals <- colSums(observed)
  grand_total <- sum(observed)
  
  expected <- outer(row_totals, col_totals) / grand_total
  
  # Calculate the chi-square statistic
  chi_square_statistic <- sum((observed - expected)^2 / expected)
  
  return(chi_square_statistic)
}

# Generate two random numeric vectors
set.seed(123)
x <- sample(1:5, 1000, replace = TRUE)
y <- sample(1:5, 1000, replace = TRUE)

# Run microbenchmark to compare chisq.test and fast_chisq_test
benchmark_results <- microbenchmark(
  chisq_test = chisq.test(x, y)$statistic,
  fast_chisq_test = fast_chisq_test(x, y),
  times = 1000
)

# Print benchmark results
print(benchmark_results)

## -----------------------------------------------------------------------------
fast_table <- function(x, y) {
  # Determine the range of possible values in x and y
  max_x <- max(x)
  max_y <- max(y)
  
  # Initialize a matrix with zeros to store counts
  counts <- matrix(0, nrow = max_x, ncol = max_y)
  
  # Populate the counts matrix
  for (i in seq_along(x)) {
    counts[x[i], y[i]] <- counts[x[i], y[i]] + 1
  }
  
  return(counts)
}

## -----------------------------------------------------------------------------
fast_chisq_test <- function(x, y) {
  # Generate the contingency table using fast_table
  observed <- fast_table(x, y)
  
  # Calculate the chi-square statistic
  row_totals <- rowSums(observed)
  col_totals <- colSums(observed)
  total <- sum(observed)
  expected <- outer(row_totals, col_totals) / total
  chisq_stat <- sum((observed - expected)^2 / expected)
  
  return(chisq_stat)
}

## -----------------------------------------------------------------------------
# Generate example data
set.seed(42)
x <- sample(1:10, 1000, replace = TRUE)
y <- sample(1:10, 1000, replace = TRUE)

# Benchmark the functions
library(microbenchmark)
microbenchmark(
  chisq_test = chisq.test(table(x, y))$statistic,
  fast_chisq_test = fast_chisq_test(x, y),
  times = 1000
)


