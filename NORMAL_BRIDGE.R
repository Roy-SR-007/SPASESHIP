# SPASESHIP: SPArSe Estimation using SHrInkage Priors

# Considering a linear regression model, y = Xbeta + eps with eps ~ N(0, sigma_sq * I)
# and the prior over beta is the Normal RIDGE prior, Bayesian RIDGE.

### Under the normal prior, the so-called “ridge” parameter was given an 
### inverse gamma prior with shape and scale parameters 10^(-3).

library(monomvn)
library(mvtnorm)
library(statmod)

##### data generation for Bayesian RIDGE and corresponding results
##### We report the model error in each case

##### Scenario 1: 5 randomly chosen components of beta are set to 1 and rest are 0, n = 50

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
  init_params = list(beta = rep(0.5, p), sigma_sq = sigma_sq)

  niter = 1000
  beta_post_est = matrix(0, ncol = p, nrow = 1)
  
  result_normal = blasso(X, y, T = niter, beta = init_params$beta, s2 = init_params$sigma_sq, case = "ridge",
                            rd = c(10^(-3), 10^(-3)), verb = 0)
  beta_post_est[1, ] = colMeans(result_normal$beta)
  beta_post_est = as.vector(beta_post_est)
  model_error_1[data] = crossprod(beta_true - beta_post_est, 
                                Sigma_X %*% as.matrix(beta_true - beta_post_est))
  
}

##### Scenario 2: 10 randomly chosen components of beta are set to 3 and rest are 0, n = 400

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
  init_params = list(beta = rep(0.5, p), sigma_sq = sigma_sq)
  
  niter = 1000
  beta_post_est = matrix(0, ncol = p, nrow = 1)
  
  result_normal = blasso(X, y, T = niter, beta = init_params$beta, s2 = init_params$sigma_sq, case = "ridge",
                         rd = c(10^(-3), 10^(-3)), verb = 0)
  beta_post_est[1, ] = colMeans(result_normal$beta)
  beta_post_est = as.vector(beta_post_est)
  model_error_2[data] = crossprod(beta_true - beta_post_est, 
                                  Sigma_X %*% as.matrix(beta_true - beta_post_est))  
}

##### Scenario 3: beta = c(0.85, 0.85, ..., 0.85), n = 400

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
  init_params = list(beta = rep(0.5, p), sigma_sq = sigma_sq)
  
  niter = 1000
  beta_post_est = matrix(0, ncol = p, nrow = 1)
  
  result_normal = blasso(X, y, T = niter, beta = init_params$beta, s2 = init_params$sigma_sq, case = "ridge",
                         rd = c(10^(-3), 10^(-3)), verb = 0)
  beta_post_est[1, ] = colMeans(result_normal$beta)
  beta_post_est = as.vector(beta_post_est)
  model_error_3[data] = crossprod(beta_true - beta_post_est, 
                                  Sigma_X %*% as.matrix(beta_true - beta_post_est))
  
}

