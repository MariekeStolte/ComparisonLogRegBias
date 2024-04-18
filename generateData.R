library(mvtnorm)
# generateData - Generates data set for logistic regression
# Idea: Draw multivariate normally distributed and generate
# - binary variables with an equal proportion of 1s and 0s by dichotomization
# of a standard normal variable at 0
# - binary variable with proportion p of 1s and proportion 1 - p of 0s by
# dichotomizing a standard normal variable at the (1 - p) quantile
# - skewed distribution by applying the exponential function (i.e. log-normal)
# Input: n - (numeric) Number of observations
#        p.nv - (numeric) Number of normally distributed covariates
#        p.bin.bal - (numeric) Number of Ber(0.5) distributed covariates
#        p.bin.unbal - (numeric) Number of Ber(prob.bin) distributed covariates
#        p.logNv - (numeric) Number of log-normally distributed covariates
#        prob.bin - (numeric) Parameter for 2nd group of binary covariates
#        cor - (numeric) pairwise correlation of the multivariate NV from which the
#                 Covariates are generated (scalar!)
#        beta - (numeric) True parameter value (-beta_0 = ... = beta_p = beta)
#        mu.nv, sd.nv - (numeric) parameters for normally distributed covariates
#        mu.logNv, sd.logNv - (numeric) parameters for log-normally distributed
#                                 covariates
# Output: data.frame with columns: y - target variable
#                                  X1 - to X(p + 1) - covariates
# Attention: sd.nv, sd.logNv and cor must fit together so that the resulting
# covariance matrix is still p.s.d!
#
#
generateData <- function(n, p.nv, p.bin.bal, p.bin.unbal, p.logNv, prob.bin, cor,
                          beta, mu.nv = 0, mu.logNv = 0, sd.nv = 1, sd.logNv = 1) {
  p <- p.nv + p.bin.bal + p.bin.unbal + p.logNv
  sigma <- matrix(cor, nrow = p, ncol = p)
  diag(sigma) <- c(rep(sd.nv^2, p.nv), rep(1, p.bin.bal + p.bin.unbal), 
                   rep(sd.logNv^2, p.logNv))
  cutoff <- qnorm(1 - prob.bin)
  X <- rmvnorm(n, mean = c(rep(mu.nv, p.nv), 
                                    rep(0, p.bin.bal + p.bin.unbal),
                                    rep(mu.logNv, p.logNv)),
                        sigma = sigma)
  if(p.bin.bal > 0) {
    X[, (p.nv + 1):(p.nv + p.bin.bal)] <- as.numeric(X[, (p.nv + 1):(p.nv + p.bin.bal)] >= 0)
  }
  if(p.bin.unbal > 0) {
    X[, (p.nv + p.bin.bal + 1):(p.nv + p.bin.bal + p.bin.unbal)] <- 
      as.numeric(X[, (p.nv + p.bin.bal + 1):(p.nv + p.bin.bal + p.bin.unbal)] >= cutoff)
  }
  if(p.logNv > 0) {
    X[, (p - p.logNv + 1):p] <- exp(X[, (p - p.logNv + 1):p])
  }
  pi <- 1 / (1 + exp(-rowSums(X * beta) + beta))
  y <- as.numeric(runif(n) <= pi)
  data <- data.frame(y = y, X)
  return(data)
}
