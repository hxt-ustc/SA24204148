## -----------------------------------------------------------------------------
set.seed(123)

# Parameters
N <- 1000  # total number of hypotheses
m <- 10000  # number of simulation replicates
alpha <- 0.1  # nominal significance level

# Function to compute FWER, FDR, and TPR for one simulation
compute_metrics <- function(p_values, true_nulls, method) {
  # Apply p-value adjustments
  if (method == "bonferroni") {
    p_adjusted <- pmin(1, p_values * N)  # Bonferroni correction
  } else if (method == "bh") {
    p_adjusted <- p.adjust(p_values, method = "BH")  # B-H correction
  }
  
  # Decision: reject if p_adjusted < alpha
  rejected <- p_adjusted < alpha
  
  # Calculate metrics
  FWER <- any(rejected[true_nulls])  # Family-Wise Error Rate (any false positive)
  FDP <- sum(rejected[true_nulls]) / max(1, sum(rejected))  # False Discovery Proportion (FDR)
  TPR <- sum(rejected[!true_nulls]) / sum(!true_nulls)  # True Positive Rate
  
  return(c(FWER = FWER, FDR = FDP, TPR = TPR))
}

# Simulation
results_bonferroni <- matrix(0, nrow = m, ncol = 3)
results_bh <- matrix(0, nrow = m, ncol = 3)

for (i in 1:m) {
  # Generate p-values
  true_nulls <- c(rep(TRUE, 950), rep(FALSE, 50))  # 950 null, 50 alternative
  p_values_null <- runif(950)  # null hypotheses: uniform distribution
  p_values_alt <- rbeta(50, 0.1, 1)  # alternative hypotheses: beta distribution
  p_values <- c(p_values_null, p_values_alt)
  
  # Bonferroni correction
  results_bonferroni[i, ] <- compute_metrics(p_values, true_nulls, method = "bonferroni")
  
  # Benjamini-Hochberg correction
  results_bh[i, ] <- compute_metrics(p_values, true_nulls, method = "bh")
}

# Calculate mean FWER, FDR, and TPR over all simulations
mean_bonferroni <- colMeans(results_bonferroni)
mean_bh <- colMeans(results_bh)

# Output the results as a 3x2 table
results_table <- matrix(c(mean_bonferroni, mean_bh), nrow = 3, byrow = FALSE)
colnames(results_table) <- c("Bonferroni correction", "B-H correction")
rownames(results_table) <- c("FWER", "FDR", "TPR")

print(results_table)

## -----------------------------------------------------------------------------
# Data: times between failures
times <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)

# MLE of lambda
lambda_mle <- 1 / mean(times)

# Bootstrap to estimate bias and standard error
set.seed(123)  # For reproducibility
n_bootstrap <- 1000  # Number of bootstrap samples
lambda_boot <- numeric(n_bootstrap)

for (i in 1:n_bootstrap) {
  boot_sample <- sample(times, replace = TRUE)
  lambda_boot[i] <- 1 / mean(boot_sample)
}

# Bias estimate
bias_estimate <- mean(lambda_boot) - lambda_mle

# Standard error estimate
se_estimate <- sd(lambda_boot)

# Results
cat("MLE of lambda:", lambda_mle, "\n")
cat("Bootstrap estimate of lambda:", mean(lambda_boot), "\n")
cat("Bootstrap bias estimate:", bias_estimate, "\n")
cat("Bootstrap standard error estimate:", se_estimate, "\n")

## -----------------------------------------------------------------------------
set.seed(123)  # For reproducibility
n <- length(times)  # Number of original observations
B <- 1000  # Number of bootstrap samples

# Store bootstrap estimates for mean time between failures (1/lambda)
mean_time_boot <- numeric(B)

for (i in 1:B) {
  boot_sample <- sample(times, n, replace = TRUE)
  lambda_boot <- 1 / mean(boot_sample)  # MLE for each bootstrap sample
  mean_time_boot[i] <- 1 / lambda_boot  # 1/lambda for each bootstrap sample
}


## -----------------------------------------------------------------------------
alpha <- 0.05
mean_time_hat <- mean(mean_time_boot)  # Mean estimate of 1/lambda
se_boot <- sd(mean_time_boot)  # Standard error of bootstrap estimates

# Standard normal confidence interval
z <- qnorm(1 - alpha / 2)
ci_normal <- c(mean_time_hat - z * se_boot, mean_time_hat + z * se_boot)

## -----------------------------------------------------------------------------
# Basic bootstrap CI
ci_basic <- 2 * mean_time_hat - quantile(mean_time_boot, c(1 - alpha / 2, alpha / 2))

## -----------------------------------------------------------------------------
# Percentile bootstrap CI
ci_percentile <- quantile(mean_time_boot, c(alpha / 2, 1 - alpha / 2))

## ----include=TRUE,echo=FALSE--------------------------------------------------
library(boot)

## -----------------------------------------------------------------------------
# Calculate BCa confidence interval for the mean time between failures
bca_ci <- boot.ci(boot(mean_time_boot, function(x, i) mean(x[i]), R = B), type = "bca")
ci_bca <- c(bca_ci$bca[4], bca_ci$bca[5])

## -----------------------------------------------------------------------------
ci_results <- data.frame(
  Method = c("Standard Normal", "Basic", "Percentile", "BCa"),
  Lower = c(ci_normal[1], ci_basic[1], ci_percentile[1], ci_bca[1]),
  Upper = c(ci_normal[2], ci_basic[2], ci_percentile[2], ci_bca[2])
)

print(ci_results)

