# SPASESHIP: SPArSe Estimation using SHrInkage Priors

# Considering a linear regression model, y = Xbeta + eps with eps ~ N(0, sigma_sq * I)
# and GDP prior over the model parameters, we execute the Gibbs sampler given in
# Armagan et al. (https://arxiv.org/pdf/1104.0861) [Section 3]

library(mvtnorm)
library(statmod)

# Griddy sampling for the hyper-parameters, alpha and eta
alpha_eta_hyper_prior_update = function(alpha_init, eta_init, beta, sigma_sq) {
  
  p = length(beta)
  a_grid = c(0.1, 0.3, 0.5, 0.7, 0.9) # grid values of a, corresponding to alpha 
  e_grid = c(0.1, 0.3, 0.5, 0.7, 0.9) # grid values of e, corresponding to eta

  w_a = c() # weights for a
  w_e = c() # weights for e
  
  for(k in 1:length(a_grid)) {
      w_a[k] = (((1-a_grid[k]) / a_grid[k]) ^ p) *
        prod((1 + abs(beta) / (sqrt(sigma_sq) * eta_init)) ^ (-(1 / a_grid[k])))
  }
  
  for(k in 1:length(e_grid)) {
    w_e[k] = ((e_grid[k] / (1 - e_grid[k])) ^ p) *
      prod((1 + (e_grid[k] * abs(beta)) / (sqrt(sigma_sq) * (1 - e_grid[k]))) ^ (-alpha_init - 1))
  }
  
  w_a_N = w_a / sum(w_a) # normalizing the weights for a
  w_e_N = w_e / sum(w_e) # normalizing the weights for e
  
  alpha = (1 / sample(a_grid, 1, prob = w_a_N)) - 1 # sample of alpha
  eta = (1 / sample(e_grid, 1, prob = w_e_N)) - 1 # sample for eta
  
  return(list(alpha = alpha, eta = eta))
}


Gibbs_Armagan = function(y, X, init_params, niter, alpha_fix = FALSE, eta_fix = FALSE) {
  
  # initialization of the parameters
  
  # regression coefficients
  beta = matrix(0, ncol = ncol(X), nrow = niter)
  beta[1, ] = init_params$beta
  # data augmentation parameters
  T_mat = matrix(0, ncol = ncol(X), nrow = niter)
  T_mat[1, ] = init_params$tau
  # model variance
  sigma_sq = c()
  sigma_sq[1] = init_params$sigma_sq
  # GDP prior: scale mixture exponential parameter 
  lambda = matrix(0, ncol = ncol(X), nrow = niter)
  lambda[1 , ] = init_params$lambda
  # GDP prior: Gamma shape parameter
  alpha = c()
  alpha[1] = init_params$alpha
  # GDP prior: Gamma scale parameter
  eta = c()
  eta[1] = init_params$eta
  
  # running the Gibbs sampler for niter iterations
  print(paste("Starting Gibbs Sampler in Armagan et al. (2001) with ", niter, " iterations"))
  for(t in 2:niter) {
    
    mean_beta = as.vector(solve(crossprod(X) + solve(diag(T_mat[t-1, ]))) %*% crossprod(X, y))
    sigma_beta = sigma_sq[t - 1] * solve(crossprod(X) + solve(diag(T_mat[t-1, ])))
    
    # update of beta according as, [beta|(sigma_sq, T, y)]
    beta[t , ] = as.vector(rmvnorm(1, mean = mean_beta, sigma = sigma_beta))
    
    shape_sigma_sq = (ncol(X) + nrow(X)) / 2
    scale_sigma_sq = (0.5 * as.numeric(crossprod(y - X %*% beta[t, ]))) + (0.5 * as.numeric(as.matrix(crossprod(beta[t, ], solve(diag(T_mat[t-1, ])))) 
                                                       %*% as.matrix(beta[t , ])))
    
    # update of sigma_sq according as, [sigma_sq|(beta, T, y)]
    sigma_sq[t] = 1 / rgamma(1, shape = shape_sigma_sq, scale = scale_sigma_sq)
    
    for(j in 1:ncol(X)) {
      
      # update of lambda's according as, [lambda_j|(beta_j, sigma_sq)]
      lambda[t, j] = rgamma(1, shape = alpha[t-1] + 1, scale = eta[t-1] + (abs(beta[t, j]) / sqrt(sigma_sq[t])))
      
      # update of tau's according as, [1 / tau_j|(beta_j, lambda_j, sigma_sq)]
      T_mat[t, j] = 1 / rinvgauss(1, mean = sqrt(((lambda[t, j] ^ 2) * sigma_sq[t]) / (beta[t , j] ^ 2)),
                                  dispersion = lambda[t, j])
      
    }
    
    if(alpha_fix == TRUE) {
      alpha[t] = 1
    }else {
      # griddy update of alpha
      alpha[t] = alpha_eta_hyper_prior_update(alpha_init = alpha[t-1],
                                              eta_init = eta[t-1], beta[t, ], sigma_sq[t])$alpha
    }
    
    if(eta_fix == TRUE) {
      eta[t] = 1
    }else {
      # griddy update of eta
      eta[t] = alpha_eta_hyper_prior_update(alpha_init = alpha[t],
                                            eta_init = eta[t-1], beta[t, ], sigma_sq[t])$eta 
    }
  }
  
  message("--------------------------------------------------------")
  
  return(list(beta = beta, sigma_sq = sigma_sq, lambda = lambda, T_mat = T_mat, alpha = alpha, eta = eta))
}


##### data generation for the above Armagan Gibbs sampler and corresponding results
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
  init_params = list(beta = rep(0.5, p), sigma_sq = sigma_sq, tau = 0.1, lambda = 1, alpha = 1, eta = 1)
  
  niter = 1000
  beta_post_est = matrix(0, ncol = p, nrow = 1)
 
  result_Gibbs_Armagan = Gibbs_Armagan(y, X, init_params, niter, alpha_fix = TRUE, eta_fix = TRUE)
  beta_post_est[1, ] = colMeans(result_Gibbs_Armagan$beta)
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
  init_params = list(beta = rep(0.5, p), sigma_sq = sigma_sq, tau = 0.1, lambda = 1, alpha = 1, eta = 1)
  
  niter = 1000
  beta_post_est = matrix(0, ncol = p, nrow = 1)
  result_Gibbs_Armagan = Gibbs_Armagan(y, X, init_params, niter, alpha_fix = TRUE, eta_fix = TRUE)
  beta_post_est[1, ] = colMeans(result_Gibbs_Armagan$beta)
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
  init_params = list(beta = rep(0.5, p), sigma_sq = sigma_sq, tau = 0.1, lambda = 1, alpha = 1, eta = 1)
  
  niter = 1000
  beta_post_est = matrix(0, ncol = p, nrow = 1)
  result_Gibbs_Armagan = Gibbs_Armagan(y, X, init_params, niter, alpha_fix = FALSE, eta_fix = FALSE)
  beta_post_est[1, ] = colMeans(result_Gibbs_Armagan$beta)
  beta_post_est = as.vector(beta_post_est)
  
  model_error_3[data] = crossprod(beta_true - beta_post_est, 
                                  Sigma_X %*% as.matrix(beta_true - beta_post_est))
  
}

