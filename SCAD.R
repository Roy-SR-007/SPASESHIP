# SPASESHIP: SPArSe Estimation using SHrInkage Priors

# Considering a linear regression model, y = Xbeta + eps with eps ~ N(0, sigma_sq * I)
# and beta is estimated using the SCAD technique. SCAD penalty.

library(SIS)
library(mvtnorm)
library(statmod)

##### data generation for SCAD penalty and corresponding results
##### We report the model error in each case.

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
  beta_true_dummy = beta_true
  X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
  #for(j in 1:ncol(X)) {
  #  X[ , j] = as.vector((X[ , j] - mean(X[ , j])) / as.numeric(sqrt(crossprod(X[ , j] - mean(X[ , j])))))
  #}
  y = as.numeric(X %*% beta_true) + error
  #y = y - mean(y)
  #init_params = list(beta = rep(0.5, p), sigma_sq = sigma_sq, tau = 0.1, lambda = 1, alpha = 1, eta = 1)
  
  niter = 1000
  result_SCAD = SIS(X, y, family = "gaussian", penalty = "SCAD", iter.max = niter)
  beta_true_dummy[result_SCAD$ix] = as.vector(result_SCAD$coef.est)[-1]
  beta_SCAD_est = beta_true_dummy
  model_error_1[data] = crossprod(beta_true - beta_SCAD_est, 
                                Sigma_X %*% as.matrix(beta_true - beta_SCAD_est))
  
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
  beta_true_dummy = beta_true
  X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
  #for(j in 1:ncol(X)) {
  #  X[ , j] = as.vector((X[ , j] - mean(X[ , j])) / as.numeric(sqrt(crossprod(X[ , j] - mean(X[ , j])))))
  #}
  y = as.numeric(X %*% beta_true) + error
  #y = y - mean(y)
  #init_params = list(beta = rep(0.5, p), sigma_sq = sigma_sq, tau = 0.1, lambda = 1, alpha = 1, eta = 1)
  
  niter = 1000
  result_SCAD = SIS(X, y, family = "gaussian", penalty = "SCAD", iter.max = niter)
  beta_true_dummy[result_SCAD$ix] = as.vector(result_SCAD$coef.est)[-1]
  beta_SCAD_est = beta_true_dummy
  model_error_2[data] = crossprod(beta_true - beta_SCAD_est, 
                                  Sigma_X %*% as.matrix(beta_true - beta_SCAD_est))
  
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
  beta_true_dummy = beta_true
  X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
  #for(j in 1:ncol(X)) {
  #  X[ , j] = as.vector((X[ , j] - mean(X[ , j])) / as.numeric(sqrt(crossprod(X[ , j] - mean(X[ , j])))))
  #}
  y = as.numeric(X %*% beta_true) + error
  #y = y - mean(y)
  #init_params = list(beta = rep(0.5, p), sigma_sq = sigma_sq, tau = 0.1, lambda = 1, alpha = 1, eta = 1)
  
  niter = 1000
  result_SCAD = SIS(X, y, family = "gaussian", penalty = "SCAD", iter.max = niter)
  beta_true_dummy[result_SCAD$ix] = as.vector(result_SCAD$coef.est)[-1]
  beta_SCAD_est = beta_true_dummy
  model_error_3[data] = crossprod(beta_true - beta_SCAD_est, 
                                  Sigma_X %*% as.matrix(beta_true - beta_SCAD_est))
  
}