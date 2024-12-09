## -----------------------------------------------------------------------------
# Function to calculate the k-th term in the series
compute_kth_term <- function(k, a, d) {
  # Euclidean norm of vector a
  norm_a <- sqrt(sum(a^2))
  
  # Calculate each component of the exponent in the term
  term <- (-1)^k * exp(
    (2 * k + 2) * log(norm_a) - 
    lgamma(k + 2) - 
    (k + 2) * log(2) + 
    lgamma((d + 1) / 2) + 
    lgamma(k + 0.5) - 
    lgamma(k + d / 2 + 1)
  )
  
  return(term)
}

## -----------------------------------------------------------------------------
compute_sum <- function(a, d, k_out = 0 ,tolerance = 1e-10, max_iter = 1000) {
  total_sum <- 0
  k <- 0
  
  repeat {
    term <- compute_kth_term(k, a, d)
    total_sum <- total_sum + term
    
    # Check for convergence
    if (abs(term) < tolerance || k >= max_iter) {
      break
    }
    
    k <- k + 1
  }
  if(k_out) cat("k=", k, "\n")
  return(total_sum)
}

## -----------------------------------------------------------------------------
compute_sum_vector <- function(a, d, k) {
  vec <- seq(k)-1
  a_norm <- sqrt(sum(a^2))
  i <- rep(c(1,-1), length = k)
  
  part1 <- lgamma(vec+1.5)-lgamma(vec+d/2+1)
  part2 <- (2*vec+2)*log(a_norm) - log(factorial(vec)) - vec*log(2) - log(2*vec+1)-log(2*vec+2)
  
  return(sum(i*exp(part1+part2))*gamma(d/2+0.5))
}

## -----------------------------------------------------------------------------
# Vector a and dimension d
a <- c(1, 2)
d <- 2

# Compute the sum
result <- compute_sum(a, d, k_out = 1)
print(paste("The sum is:", result))

compute_sum_vector(a, d, 17)

## -----------------------------------------------------------------------------
d = 10
a = seq(d)

system.time(
  {
    for (i in seq(1000)) {
      compute_sum(a, d)
    }
  }
)

system.time(
  {
    for (i in seq(1000)) {
      compute_sum_vector(a, d, 17)
    }
  }
)

## -----------------------------------------------------------------------------
# Define a function to calculate Sk-1(a) - Sk(a)
compute_difference <- function(a, k) {
  # Calculate Sk-1(a)
  critical_value_k_minus_1 <- sqrt(a^2 * (k - 1) / (k - a^2))
  prob_k_minus_1 <- pt(critical_value_k_minus_1, df = k - 1, lower.tail = FALSE)
  
  # Calculate Sk(a)
  critical_value_k <- sqrt(a^2 * k / (k + 1 - a^2))
  prob_k <- pt(critical_value_k, df = k, lower.tail = FALSE)
  
  # Return the difference
  return(prob_k_minus_1 - prob_k)
}

# Set the value of k, for example, k = 4
k_set <- c(4 : 25, 100, 500, 1000)

result1 <- NULL

for (k in k_set) {
  # Use uniroot to find the value of a in the interval (0, sqrt(k))
  result <- uniroot(compute_difference, interval = c(1e-5, ifelse(k>100,sqrt(k)-4,sqrt(k)-1e-5)), k = k)
  a_intersection <- result$root
  result1 <- c(result1, a_intersection)

  cat("Intersection point A(k) for k =", k, "is a =", a_intersection, "\n")
}

## -----------------------------------------------------------------------------
# Define the function to compute the integrals and find 'a'
solve_for_a <- function(k) {
  # Define the target function based on the equation
  target_function <- function(a) {
    # Compute ck and ck-1
    ck <- sqrt(a^2 * k / (k + 1 - a^2))
    ck_minus_1 <- sqrt(a^2 * (k - 1) / (k - a^2))
    
    # Compute the left-hand side integral
    left_integral <- integrate(function(u) (1 + u^2 / (k - 1))^(-k / 2), 0, ck_minus_1)$value
    left_side <- 2 * exp(lgamma(k / 2) - lgamma((k - 1) / 2)) / (sqrt( (k - 1))) * left_integral
    
    # Compute the right-hand side integral
    right_integral <- integrate(function(u) (1 + u^2 / k)^(-(k + 1) / 2), 0, ck)$value
    right_side <- 2 * exp(lgamma((k + 1) / 2)-lgamma(k / 2)) / (sqrt( k) )  * right_integral
    
    # Return the difference
    return(left_side - right_side)
  }
  
  # Use a root-finding method to find 'a' such that target_function(a) = 0
  result <- uniroot(target_function, c(1, 1.8))$root
  return(result)
  
}

result2 <- NULL

# Example usage with a specific value of k
k_set <- c(4 : 25, 100, 500, 1000)
for (k in k_set) {
  a_solution <- solve_for_a(k)
  result2 <- c(result2, a_solution)
  cat("the root for k =",k , "is", a_solution,"\n")
}

## -----------------------------------------------------------------------------
result1-result2

## -----------------------------------------------------------------------------
y = c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau = 1
lambda1 = sum(y<tau)/sum(y)
cat("The estimate of lambda by mle is", lambda1, "\n")

## -----------------------------------------------------------------------------
# E-M algorithm
tol <- 1e-6  # Convergence tolerance
max_iter <- 1000  # Maximum number of iterations
lambda = lambda1
for (iter in 1:max_iter) {
  # E-step: Calculate expected T_i for censored and uncensored values
  expected_T <- ifelse(y < tau, y, tau + 1 / lambda)
  
  # M-step: Update lambda
  new_lambda <- length(y) / sum(expected_T)
  
  # Check convergence
  if (abs(new_lambda - lambda) < tol) {
    lambda <- new_lambda
    break
  }
  
  # Update lambda for the next iteration
  lambda <- new_lambda
}

# Output the estimated lambda
lambda

