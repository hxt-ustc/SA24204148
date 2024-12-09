## ----include=TRUE, echo=FALSE, warning=FALSE, message=FALSE-------------------
library(bootstrap)
library(GGally)
library(boot)
library(DAAG)
library(ggplot2)

## -----------------------------------------------------------------------------
# The dataset contains the test scores in five subjects
data(scor)

head(scor)

# Panel display of scatter plots with correlations
ggpairs(scor, title = "Scatter Plot Matrix with Correlations")


## -----------------------------------------------------------------------------
# Function to compute the proportion of variance explained by the first principal component
pve_first_pc <- function(data, indices, loc = NULL) {
  d <- data[indices, ]  
  cov_matrix <- cov(d)
  eigenvalues <- eigen(cov_matrix)$values
  theta_hat <- eigenvalues[1] / sum(eigenvalues)
  return(theta_hat)
}

# Number of observations
n <- nrow(scor)

# Placeholder for jackknife estimates of the correlation
jackknife_theta <- numeric(n)

# Calculate correlation leaving one observation out at a time
for (i in 1:n) {
  # Leave the i-th observation out and calculate the correlation
  jackknife_theta[i] <- pve_first_pc(scor,(seq(n)[-i]))
}

# Jackknife estimate of the correlation (mean of leave-one-out correlations)
jackknife_mean <- mean(jackknife_theta)

# original theta_hat

theta_full <- pve_first_pc(scor,seq(n))

# Jackknife bias estimate
bias_jack <- (n - 1) * (jackknife_mean - theta_full)

# Jackknife standard error estimate
se_jack <- sqrt(((n - 1) / n) * sum((jackknife_theta - jackknife_mean)^2))

# Output the results
cat("Original Theta (Full Data):", theta_full, "\n")
cat("Jackknife Estimate:", jackknife_mean, "\n")
cat("Jackknife Bias Estimate:", bias_jack, "\n")
cat("Jackknife Standard Error Estimate:", se_jack, "\n")

## -----------------------------------------------------------------------------
# Bootstrap to estimate bias and standard error for proportion of variance explained
set.seed(123)
boot_pve_results <- boot(data = scor, statistic = pve_first_pc, R = 1000, loc = NULL)

# Standard error and bias
boot_pve_bias <- mean(boot_pve_results$t) - boot_pve_results$t0
boot_pve_se <- sd(boot_pve_results$t)

result <- rbind(c(boot_pve_bias,boot_pve_se),c(bias_jack,se_jack))

colnames(result) <- c("bias", "se")
rownames(result) <- c("bootstrap", "jackknife")

round(result,4)

## -----------------------------------------------------------------------------
data(ironslag)
# par(mfrow=c(2,2), mar = c(5, 4, 4, 2) + 0.1)

a <- seq(10, 40, .1) #sequence for plotting fits
L1 <- lm(magnetic ~ chemical, data = ironslag)
plot(ironslag$chemical, ironslag$magnetic, main = "Linear", pch = 16 ,xlab = "")
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)

L2 <- lm(magnetic ~ chemical + I(chemical^2), data = ironslag)
plot(ironslag$chemical, ironslag$magnetic, main="Quadratic", pch=16, xlab = "", ylab = "")
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)

L3 <- lm(log(magnetic) ~ chemical, data = ironslag)
plot(ironslag$chemical, ironslag$magnetic, main="Exponential", pch=16,)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)

L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3), data = ironslag)
plot(ironslag$chemical, ironslag$magnetic, main="Cubic", pch=16, ylab = "")
yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2 + L4$coef[4] * a^3
lines(a, yhat4, lwd=2)

## -----------------------------------------------------------------------------
n <- length(ironslag$magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)

# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- ironslag$magnetic[-k]
  x <- ironslag$chemical[-k]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * ironslag$chemical[k]
  e1[k] <- ironslag$magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * ironslag$chemical[k] + J2$coef[3] * ironslag$chemical[k]^2
  e2[k] <- ironslag$magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * ironslag$chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- ironslag$magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * ironslag$chemical[k] + J4$coef[3] * ironslag$chemical[k]^2 + J4$coef[4] * ironslag$chemical[k]^3
  e4[k] <- ironslag$magnetic[k] - yhat4
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
summary(L1)$adj.r.squared
summary(L2)$adj.r.squared
summary(L3)$adj.r.squared
summary(L4)$adj.r.squared

## -----------------------------------------------------------------------------
data("chickwts")
head(chickwts)
boxplot(formula(chickwts),data = chickwts)

## -----------------------------------------------------------------------------
x <- sort(as.vector(chickwts[chickwts$feed == "soybean",1]))
y <- sort(as.vector(chickwts[chickwts$feed == "linseed",1]))

# Define the Cramer-von Mises statistic based on the provided formula
cvm_statistic <- function(x, y) {
  # Sample sizes
  n <- length(x)
  m <- length(y)
  
  # Combine x and y into a single vector
  combined <- c(x, y)
  
  # Sort combined data
  combined_sorted <- sort(combined)
  
  # ECDF for x and y
  Fn <- ecdf(x)  # ECDF for sample x
  Gm <- ecdf(y)  # ECDF for sample y
  
  # Cramer-von Mises test statistic as per the given formula
  W2 <- (n * m) / (n + m)^2 * (
    sum((Fn(x) - Gm(x))^2) +  # Summation for x_i terms
    sum((Fn(y) - Gm(y))^2)    # Summation for y_j terms
  )
  
  return(W2)
}

# Perform the permutation test
permutation_test <- function(x, y, n_perm = 1000) {
  # Observed Cramer-von Mises statistic
  observed_stat <- cvm_statistic(x, y)
  
  # Combine x and y
  combined <- c(x, y)
  n <- length(combined)
  
  # Initialize a vector to store permutation statistics
  perm_stats <- numeric(n_perm)
  
  # Run permutation
  for (i in 1:n_perm) {
    # Permute the combined data
    permuted <- sample(combined, size = n, replace = FALSE)
    
    # Split permuted data into two groups
    perm_x <- permuted[1:length(x)]
    perm_y <- permuted[(length(x) + 1):n]
    
    # Calculate the Cramer-von Mises statistic for the permuted samples
    perm_stats[i] <- cvm_statistic(perm_x, perm_y)
  }
  
  # Calculate p-value
  p_value <- mean(perm_stats >= observed_stat)
  
  # Return observed statistic and p-value
  return(list(observed_stat = observed_stat, p_value = p_value, perm_stats = perm_stats))
}

# Example usage
set.seed(123)

# Perform the permutation test
result <- permutation_test(x, y, n_perm = 1000)

# Display the results
result$observed_stat  # Observed Cramer-von Mises statistic
result$p_value        # Permutation test p-value

## -----------------------------------------------------------------------------
hist(result$perm_stats, main = "", xlab = "perm_stats(p=0.422)",freq = FALSE, breaks = "scott")
points(result$observed_stat, 0, cex = 1, pch = 16) #observed W2

## -----------------------------------------------------------------------------
set.seed(123)  # For reproducibility

# Generate two random samples
x <- rnorm(30)
y <- rnorm(30)
plot(x,y)

## -----------------------------------------------------------------------------
# 1. Compute the observed Spearman correlation
observed_corr <- cor(x, y, method = "spearman")

# 2. Perform the permutation test
perm_test_spearman <- function(x, y, num_permutations = 1000) {
  # Initialize a vector to store permuted correlations
  permuted_corrs <- numeric(num_permutations)
  combined <- c(x, y)
  
  n = length(x)
  m = length(y)
  
  # Loop over the number of permutations
  for (i in 1:num_permutations) {
    
    permuted_index <- sample(1:(n+m), n, FALSE)
    
    # Calculate the Spearman correlation with shuffled data
    permuted_corrs[i] <- cor(combined[permuted_index], combined[-permuted_index], method = "spearman")
  }
  
  # Return the permuted correlations
  return(permuted_corrs)
}

# Number of permutations
num_permutations <- 1000

# Run the permutation test
permuted_corrs <- perm_test_spearman(x, y, num_permutations)

# 3. Calculate the permutation p-value
p_value_perm <- mean(abs(permuted_corrs) >= abs(observed_corr))

# 4. Compare with the p-value from cor.test
spearman_test <- cor.test(x, y, method = "spearman")

# Results
cat("Observed Spearman correlation:", observed_corr, "\n")
cat("Permutation p-value:", p_value_perm, "\n")
cat("Spearman test p-value:", spearman_test$p.value, "\n")

# Plot histogram of permuted correlations
hist(permuted_corrs, main = "Permutation Test: Spearman Correlations",
     xlab = "Spearman Correlation", col = "lightblue", breaks = 30)
abline(v = observed_corr, col = "red", lwd = 2, lty = 2)  # Mark observed correlation

