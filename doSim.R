source("repSim.R")
# doSim - Performs simulation for desired parameter settings
# Input: ns - (numeric) vector with settings for n
#        betas - (numeric) vector with settings for beta
#        cors - (numeric) vector with settings for correlation of the
#                 Covariates
#        scenarios - (data.frame) with columns p, type (name of scenario),
#                       p.nv (number N(0,1) covariates), p.bin.bal (number
#                       Ber(0.5) covariates), p.bin.unbal (number Ber(0.2)
#                       covariates), p.logNv (number logNv(-1, 0.75)
#                       Covariates) -> Contains settings for all
#                       distributions of the design matrix to be considered
#        repl - (numeric) number of repetitions (default: 1000)
#        seed - (numeric) Seed (default: 24092021)
# Output: data.frame with columns estimate, p.value (p-value for
#           Test for beta_i = 0), beta (True parameter value), epv (Events per
#           Variable for underlying data set), angle (angle of
#           true parameter vector to the EV of the largest EW of X^TX),
# separation (Has separability been determined?), estimator
# (estimator used), n, p, cor, scenario (see input), distr
# (Distribution of the associated covariable or intercept, if
# Intercept), replication (repetition for exactly these parameters-
# settings), ID.comb (consecutive number, one for each parameter
# combination), correlation.1,...correlation.Xp observed
# Correlations
#
# calcParamsLogNv <- function(mu.logNv, sd.logNv) {
#   s <- sqrt(log((sd.logNv / mu.logNv)^2 + 1))
#   m <- log(mu.logNv) - s^2 / 2
#   return(c(m, s))
# }
# params <- calcParamsLogNv(1, 1)

doSim <- function(ns, betas, cors, scenarios, p = 1, repl1 = 1000, repl2 = 5000,
                  seed = 24092021, prob.bin = 0.2, 
                  mu.logNv = (-log(exp(1) - 1) - 1) / 2, sd.logNv = 1) {
  set.seed(seed)
  scenarios <- scenarios[scenarios$p == p, ]
  combs <- expand.grid(n = ns, beta = betas, cor = cors, scen.ind = 1:nrow(scenarios))
  # For p = 1 no different correlations -> delete corresponding rows:
  combs <- combs[!(p == 1 & combs$cor != cors[1]), ]
  res <- mapply(repSim, n = combs$n, beta = combs$beta, cor = combs$cor,
                ind.scen = combs$scen.ind, 
                MoreArgs = list(repl1 = repl1, repl2 = repl2, scenarios = scenarios, 
                                mu.logNv = mu.logNv, sd.logNv = sd.logNv, 
                                prob.bin = prob.bin), 
                SIMPLIFY = FALSE)
  ls <- sapply(res, nrow)
  res <- do.call(rbind, res)
  res$ID.comb <- rep(1:nrow(combs), times = ls)
  return(res)
}
