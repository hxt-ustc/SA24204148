## ----include=FALSE------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

## -----------------------------------------------------------------------------
cppFunction('
DataFrame gibbs_sampler_rcpp(int iterations, int n, double a, double b, double x_init, double y_init) {
  // Storage for samples
  NumericVector x_samples(iterations);
  NumericVector y_samples(iterations);
  
  // Initialize x and y
  double x = x_init;
  double y = y_init;
  
  // Gibbs Sampling
  for (int i = 0; i < iterations; i++) {
    // Step 1: Sample x from Binomial(n, y)
    x = R::rbinom(n, y);
    
    // Step 2: Sample y from Beta(x + a, n - x + b)
    y = R::rbeta(x + a, n - x + b);
    
    // Store samples
    x_samples[i] = x;
    y_samples[i] = y;
  }
  
  // Return samples as a DataFrame
  return DataFrame::create(
    Named("x") = x_samples,
    Named("y") = y_samples
  );
}
')

## -----------------------------------------------------------------------------
gibbs_sampler_r <- function(iterations, n, a, b, x_init, y_init) {
  # Storage for samples
  x_samples <- numeric(iterations)
  y_samples <- numeric(iterations)
  
  # Initialize x and y
  x <- x_init
  y <- y_init
  
  # Gibbs Sampling
  for (i in 1:iterations) {
    # Step 1: Sample x from Binomial(n, y)
    x <- rbinom(1, size = n, prob = y)
    
    # Step 2: Sample y from Beta(x + a, n - x + b)
    y <- rbeta(1, shape1 = x + a, shape2 = n - x + b)
    
    # Store samples
    x_samples[i] <- x
    y_samples[i] <- y
  }
  
  # Return samples as a data frame
  data.frame(x = x_samples, y = y_samples)
}

## -----------------------------------------------------------------------------
# Generate samples using both methods
set.seed(123)
iterations <- 10000
n <- 10
a <- 2
b <- 2
x_init <- 5
y_init <- 0.5

samples_r <- gibbs_sampler_r(iterations, n, a, b, x_init, y_init)
samples_rcpp <- gibbs_sampler_rcpp(iterations, n, a, b, x_init, y_init)

op <- par(mfrow = c(1, 2))

# Q-Q Plot for comparing `y` samples
qqplot(
  samples_r$y, samples_rcpp$y,
  main = "Q-Q Plot: R vs. Rcpp (y samples)",
  xlab = "Quantiles of R Implementation",
  ylab = "Quantiles of Rcpp Implementation",
  col = "blue", pch = 16
)
abline(0, 1, col = "red", lwd = 2)  # Add a 45-degree line for reference

# Q-Q Plot for comparing `x` samples
qqplot(
  samples_r$x, samples_rcpp$x,
  main = "Q-Q Plot: R vs. Rcpp (x samples)",
  xlab = "Quantiles of R Implementation",
  ylab = "Quantiles of Rcpp Implementation",
  col = "blue", pch = 16
)
abline(0, 1, col = "red", lwd = 2)


## -----------------------------------------------------------------------------
# Set parameters for the Gibbs sampler
iterations <- 10000
n <- 10
a <- 2
b <- 2
x_init <- 5
y_init <- 0.5

# Benchmark the two functions
benchmark_results <- microbenchmark(
  R_implementation = gibbs_sampler_r(iterations, n, a, b, x_init, y_init),
  Rcpp_implementation = gibbs_sampler_rcpp(iterations, n, a, b, x_init, y_init),
  times = 10  # Number of repetitions
)

# Print benchmark results
print(benchmark_results)

## -----------------------------------------------------------------------------
# 运行采样器一次
samples_r <- gibbs_sampler_r(iterations, n, a, b, x_init, y_init)
samples_rcpp <- gibbs_sampler_rcpp(iterations, n, a, b, x_init, y_init)

# 绘制散点图
par(mfrow = c(1, 2))
plot(samples_r$x, samples_r$y, col = rgb(0, 0, 1, 0.3), pch = 16,
     main = "纯 R 实现: x 与 y 联合分布", xlab = "x", ylab = "y")
plot(samples_rcpp$x, samples_rcpp$y, col = rgb(1, 0, 0, 0.3), pch = 16,
     main = "Rcpp 实现: x 与 y 联合分布", xlab = "x", ylab = "y")
par(mfrow = c(1, 1))

