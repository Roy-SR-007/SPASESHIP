# Helper function to compute the stable inverse of a matrix using the Cholesky decomposition
stablesolve = function(x) {chol2inv(chol(x))}

# Helper function to compute the Euclidean norm of a vector
enorm = function(x) {return(sqrt(sum(x^2)))}

# TAVIE for location scale families: GDP Prior
#
# derivlogmother is a function that returns the derivative of the log of the base (mother) distribution of the 
# location scale family. Assume it can only evaluate one value at a time.
TAVIE_GDP = function(X, y, alpha, eta,
                     tol = 1e-5)
{
  n = length(y)
  p = ncol(X)
  xi_prev = rep(0, p)
  xi = rep(1, p)
  
  XtX = crossprod(X)
  Xty = as.numeric(crossprod(X, y))
  
  iter = 0
  while(enorm(xi - xi_prev) > tol)
  {
    Sigma = stablesolve(diag(1 / xi) + XtX)
    mu = as.numeric(Sigma %*% Xty)
    K = n / (as.numeric(sum(y ^ 2) - crossprod(Xty, Sigma %*% Xty)))
    xi_prev = xi
    kappa = as.numeric(diag(Sigma) + ((mu ^ 2) * K))
    xi = (sqrt(kappa) * (eta + sqrt(kappa))) / (alpha + 1)
    iter = iter + 1
  }
  message(paste0("Conveged in ", iter, " iterations. "))
  #message(paste0("sigma^2 = ",  K/(a+n-1)))
  return(list(xi = xi, mu = mu, Sigma = Sigma, K = K))
}


##### data generation for the above Armagan Gibbs sampler and corresponding results
##### We report the model error in each case

##### Scenario 1: 5 randomly chosen components of beta are set to 1 and rest are 0, n = 50

library(mvtnorm)
n = 50
p = 20
sigma_sq = 3 ^ 2

seed = 633
set.seed(seed)

# x_{ij}'s are standard normal with Cov(x_j, x_{j'}) = 0.5 ^ (|j - j'|)
Sigma_X = matrix(0, ncol = p, nrow = p)
for(i in 1:p) {
  for(j in 1:p) {
    Sigma_X[i, j] = 0.5 ^ (abs(i - j))
  }
}

model_error_1 = c()

# simulations for 10 data sets
for(data in 1:100) {
  
  error = rnorm(n, mean = 0, sd = sqrt(sigma_sq))
  beta_true = rep(0, p)
  beta_true[sample(1:p, 5)] = 1
  X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
  #for(j in 1:ncol(X)) {
  #  X[ , j] = as.vector((X[ , j] - mean(X[ , j])) / as.numeric(sqrt(crossprod(X[ , j] - mean(X[ , j])))))
  #}
  y = as.numeric(X %*% beta_true) + error
  #y = y - mean(y)
  
  #niter = 1000
  beta_post_est = matrix(0, ncol = p, nrow = 1)
  
  result_GDP = TAVIE_GDP(X, y, alpha = 1, eta = 1)
  
  beta_post_est[1, ] = result_GDP$mu
  beta_post_est = as.vector(beta_post_est)
  model_error_1[data] = crossprod(beta_true - beta_post_est, 
                                  Sigma_X %*% as.matrix(beta_true - beta_post_est))
  
}

##### Scenario 2: 10 randomly chosen components of beta are set to 3 and rest are 0, n = 400

library(mvtnorm)
n = 400
p = 20
sigma_sq = 3 ^ 2

seed = 633
set.seed(seed)

# x_{ij}'s are standard normal with Cov(x_j, x_{j'}) = 0.5 ^ (|j - j'|)
Sigma_X = matrix(0, ncol = p, nrow = p)
for(i in 1:p) {
  for(j in 1:p) {
    Sigma_X[i, j] = 0.5 ^ (abs(i - j))
  }
}

model_error_2 = c()

# simulations for 10 data sets
for(data in 1:100) {
  
  error = rnorm(n, mean = 0, sd = sqrt(sigma_sq))
  beta_true = rep(0, p)
  beta_true[sample(1:p, 10)] = 3
  X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
  #for(j in 1:ncol(X)) {
  #  X[ , j] = as.vector((X[ , j] - mean(X[ , j])) / as.numeric(sqrt(crossprod(X[ , j] - mean(X[ , j])))))
  #}
  y = as.numeric(X %*% beta_true) + error
  #y = y - mean(y)
  
  #niter = 1000
  beta_post_est = matrix(0, ncol = p, nrow = 1)
  
  result_GDP = TAVIE_GDP(X, y, alpha = 1, eta = 1)
  
  beta_post_est[1, ] = result_GDP$mu
  beta_post_est = as.vector(beta_post_est)
  model_error_2[data] = crossprod(beta_true - beta_post_est, 
                                  Sigma_X %*% as.matrix(beta_true - beta_post_est))
  
}

##### Scenario 3: beta = c(0.85, 0.85, ..., 0.85), n = 400

library(mvtnorm)
n = 400
p = 20
sigma_sq = 3 ^ 2

seed = 633
set.seed(seed)

# x_{ij}'s are standard normal with Cov(x_j, x_{j'}) = 0.5 ^ (|j - j'|)
Sigma_X = matrix(0, ncol = p, nrow = p)
for(i in 1:p) {
  for(j in 1:p) {
    Sigma_X[i, j] = 0.5 ^ (abs(i - j))
  }
}

model_error_3 = c()

# simulations for 10 data sets
for(data in 1:100) {
  
  error = rnorm(n, mean = 0, sd = sqrt(sigma_sq))
  beta_true = rep(0.85, p)
  X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
  #for(j in 1:ncol(X)) {
  #  X[ , j] = as.vector((X[ , j] - mean(X[ , j])) / as.numeric(sqrt(crossprod(X[ , j] - mean(X[ , j])))))
  #}
  y = as.numeric(X %*% beta_true) + error
  #y = y - mean(y)
  
  #niter = 1000
  beta_post_est = matrix(0, ncol = p, nrow = 1)
  
  result_GDP = TAVIE_GDP(X, y, alpha = 1, eta = 1)
  
  beta_post_est[1, ] = result_GDP$mu
  beta_post_est = as.vector(beta_post_est)
  model_error_3[data] = crossprod(beta_true - beta_post_est, 
                                  Sigma_X %*% as.matrix(beta_true - beta_post_est))
  
}
