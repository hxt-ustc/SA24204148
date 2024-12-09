## ----include=TRUE, echo=FALSE, warning=FALSE----------------------------------
library(coda)
library(ggplot2)
library(gridExtra)

## -----------------------------------------------------------------------------
# Set parameters
n_chains <- 3        # Number of chains
n_samples <- 5000    # Samples per chain
burn_in <- 1000      # Burn-in period
theta <- 1           # Scale parameter for standard Cauchy
eta <- 0             # Location parameter for standard Cauchy

# Metropolis-Hastings function to sample from standard Cauchy
mh_cauchy <- function(n, theta = 1, eta = 0) {
  samples <- numeric(n)
  samples[1] <- 0  # Initial value
  for (i in 2:n) {
    proposal <- samples[i - 1] + rnorm(1, mean = 0, sd = 1)
    acceptance_ratio <- dcauchy(proposal, location = eta, scale = theta) /
                        dcauchy(samples[i - 1], location = eta, scale = theta)
    if (runif(1) < acceptance_ratio) {
      samples[i] <- proposal
    } else {
      samples[i] <- samples[i - 1]
    }
  }
  return(samples)
}

## -----------------------------------------------------------------------------
# Function to calculate Gelman-Rubin diagnostic
gelman_rubin <- function(chains) {
  M <- length(chains)                   # Number of chains
  N <- length(chains[[1]])               # Length of each chain
  
  # Calculate mean and variance within each chain
  chain_means <- sapply(chains, mean)
  chain_variances <- sapply(chains, var)/N*(N+1)
  
  # Calculate the overall mean of all chains
  overall_mean <- mean(unlist(chains))
  
  # Calculate B (between-chain variance)
  B <- N * var(chain_means)
  
  # Calculate W (within-chain variance)
  W <- mean(chain_variances)
  
  # Calculate V_hat (estimated variance)
  V_hat <- (N - 1) / N * W + B / N
  
  # Calculate R_hat
  R_hat <- V_hat / W
  return(R_hat)
}

## -----------------------------------------------------------------------------
# Function to run multiple chains and check convergence
run_chains <- function(n_chains, n_samples, burn_in) {
  # Run each chain independently
  chains <- lapply(1:n_chains, function(x) mh_cauchy(n_samples))
  
  # Discard burn-in samples
  chains <- lapply(chains, function(chain) chain[(burn_in + 1):n_samples])
  
  # Calculate Gelman-Rubin diagnostic
  r_hat <- gelman_rubin(chains)
  
  return(list(r_hat = r_hat, chains = chains))
}

set.seed(12345)

# Run chains and monitor convergence
converged <- FALSE
iteration <- 1
while (!converged) {
  result <- run_chains(n_chains, n_samples, burn_in)
  r_hat <- result$r_hat
  
  cat("Iteration:", iteration, "R-hat:", r_hat, "\n")
  
  # Check if R-hat is below the threshold
  if (r_hat < 1.2) {
    converged <- TRUE
    cat("Chains have converged with R-hat <", r_hat, "\n")
  } else {
    iteration <- iteration + 1
    n_samples <- n_samples + 1000  # Increase sample size if not converged
  }
}

# Analyze the resulting chains
mcmc_samples <- unlist(result$chains)
generated_deciles <- quantile(mcmc_samples, probs = seq(0.1, 0.9, by = 0.1))
theoretical_deciles <- qcauchy(seq(0.1, 0.9, by = 0.1))

# Print results
print("Generated Deciles:")
print(generated_deciles)
print("Theoretical Deciles:")
print(theoretical_deciles)

par(mfrow = c(1, 2))
plot(mcmc_samples, type="l")
hist(mcmc_samples, breaks="scott", main="", xlab="", freq=FALSE, col = "white")
curve(dcauchy, add = T)

## -----------------------------------------------------------------------------
set.seed(123)

# Parameters
a <- 2     # Parameter a
b <- 3     # Parameter b
n <- 10    # Parameter n
iterations_per_batch <- 1000  # Number of iterations to generate per batch for each chain
chains <- 3          # Number of chains
converged <- FALSE   # Convergence flag

# Initialize storage for each chain's samples
x_samples <- list()
y_samples <- list()
for (i in 1:chains) {
  x_samples[[i]] <- numeric()
  y_samples[[i]] <- numeric()
}

# Gelman-Rubin calculation function
gelman_rubin <- function(chains_matrix) {
  m <- ncol(chains_matrix)  # Number of chains
  n <- nrow(chains_matrix)  # Number of samples per chain
  
  # Calculate mean for each chain and overall mean
  chain_means <- colMeans(chains_matrix)
  overall_mean <- mean(chains_matrix)
  
  # Calculate B (between-chain variance) and W (within-chain variance)
  B <- n * var(chain_means)
  W <- mean(apply(chains_matrix, 2, var))
  
  # Estimate of marginal posterior variance
  var_hat <- ((n - 1) / n) * W + (1 / n) * B
  
  # Gelman-Rubin statistic
  R_hat <- var_hat / W
  return(R_hat)
}

# Gibbs sampling loop until convergence
while (!converged) {
  for (chain in 1:chains) {
    # Set different initial values for each chain if starting a new run
    if (length(x_samples[[chain]]) == 0) {
      x <- sample(0:n, 1)
      y <- runif(1)
    } else {
      x <- x_samples[[chain]][length(x_samples[[chain]])]
      y <- y_samples[[chain]][length(y_samples[[chain]])]
    }
    
    # Generate samples for the current batch
    for (i in 1:iterations_per_batch) {
      # Conditional sample x | y
      x <- rbinom(1, n, y)
      
      # Conditional sample y | x
      y <- rbeta(1, x + a, n - x + b)
      
      # Store samples
      x_samples[[chain]] <- c(x_samples[[chain]], x)
      y_samples[[chain]] <- c(y_samples[[chain]], y)
    }
  }
  
  # Convert lists to matrices for Gelman-Rubin calculation
  x_matrix <- do.call(cbind, lapply(x_samples, function(x) tail(x, iterations_per_batch)))
  y_matrix <- do.call(cbind, lapply(y_samples, function(y) tail(y, iterations_per_batch)))
  
  # Calculate Gelman-Rubin statistics
  R_hat_x <- gelman_rubin(x_matrix)
  R_hat_y <- gelman_rubin(y_matrix)
  
  cat("Current R-hat for x:", R_hat_x, "\n")
  cat("Current R-hat for y:", R_hat_y, "\n")
  
  # Check convergence
  if (R_hat_x < 1.2 && R_hat_y < 1.2) {
    converged <- TRUE
    cat("Chains have converged with R-hat < 1.2\n")
  }
}

## -----------------------------------------------------------------------------
# Prepare data for plotting
# Flatten and store all samples along with their chain labels and iterations
x_df <- do.call(rbind, lapply(1:chains, function(chain) {
  data.frame(iter = seq_along(x_samples[[chain]]), 
             chain = factor(chain), 
             x = x_samples[[chain]])
}))

y_df <- do.call(rbind, lapply(1:chains, function(chain) {
  data.frame(iter = seq_along(y_samples[[chain]]), 
             chain = factor(chain), 
             y = y_samples[[chain]])
}))

p1 <- ggplot(x_df, aes(x = iter, y = x, color = chain)) +
  geom_line() + theme_minimal() + labs(title = "Trace Plot for x", x = "Iteration", y = "x") +
  scale_color_discrete(name = "Chain")

p2 <- ggplot(y_df, aes(x = iter, y = y, color = chain)) +
  geom_line() + theme_minimal() + labs(title = "Trace Plot for y", x = "Iteration", y = "y") +
  scale_color_discrete(name = "Chain")

# Plot histograms of x and y samples from all chains
p3 <- ggplot(x_df, aes(x = x, fill = chain)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  theme_minimal() + labs(title = "Histogram of x", x = "x", y = "Frequency") +
  scale_fill_discrete(name = "Chain")

p4 <- ggplot(y_df, aes(x = y, fill = chain)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  theme_minimal() + labs(title = "Histogram of y", x = "y", y = "Frequency") +
  scale_fill_discrete(name = "Chain")

# Display all plots in a grid
grid.arrange(p1, p2, p3, p4, ncol = 2)


