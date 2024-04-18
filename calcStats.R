# calcStats - Calculates target values from simulation results
# Input: res.sim - (data.frame) Result from doSim
# Output: (data.frame) Columns:
#           beta, estimator, n, p, cor, scenario, distr - As in input
#           epv, angle, separation - input averaged per parameter setting (ID.comb)
#           a.bias - Average absolute bias for this parameter setting
#           rel.bias - Relative bias (NA for beta = 0)
#           mse - Mean Squared Error
#           Power - Proportion of rejected tests on beta_i = 0
#
calcStats <- function(res.sim) {
  res.sim <- split(res.sim, res.sim$ID.comb)
  stats <- lapply(res.sim, function(dat) {
    sep <- mean(dat$separation, na.rm = TRUE)
    dat$converged[dat$estimator == "WA"] <- dat$converged[dat$estimator == "ML"]
    conv <- aggregate(dat$converged, by = list(dat$estimator), mean)
    converged <- conv$x
    names(converged) <- conv$Group.1
    converged <- converged[-which(names(converged) == "Separation")]
    tmp <- dat[dat$replication == 1 & dat$estimator != "Separation", 
               c("beta", "estimator", "n", "p", "cor",  "scenario", "distr")]
    dat <- dat[!dat$separation & !is.na(dat$separation) & 
                 dat$estimator != "Separation", ]
    dat$estimate[!(dat$converged | is.na(dat$converged)) | abs(dat$estimate) > 10] <- NA
    if(nrow(dat) == 0) {
      cor <- data.frame(t(rep(NA, sum(grepl("correlation", colnames(dat))))))
      colnames(cor) <- colnames(dat[, grepl("correlation", colnames(dat))])
      return(data.frame(tmp, angle = NA, epv = NA, bias = NA, sd = NA, 
                        median = NA, min = NA, max = NA, a.bias = NA, 
                        rel.bias = NA, mse = NA, mcse.mse = NA, power = NA, 
                        mcse.power = NA, sep = sep, converged = NA, 
                        cor))
    }
    reps <- length(unique(dat$replication))
    angle <- mean(dat$angle)
    epv <- mean(dat$epv, na.rm = TRUE)
    # 1 line corresponds to 1 repetition:
    abw <- matrix(dat$estimate - dat$beta, nrow = reps, byrow = TRUE)
    bias <- colMeans(abw, na.rm = TRUE)
    a.bias <- abs(bias)
    if(any(tmp$beta != 0)) {
      rel.bias <- a.bias / abs(tmp$beta)
    } else {
      rel.bias <- NA
    }
    mse <- colMeans(abw^2, na.rm = TRUE)
    sd <- apply(abw, 2, sd, na.rm = TRUE)
    median <- apply(abw, 2, median, na.rm = TRUE)
    min <- apply(abw, 2, min, na.rm = TRUE)
    max <- apply(abw, 2, max, na.rm = TRUE)
    power <- colMeans(matrix(dat$p.value < 0.05, nrow = reps, byrow = TRUE),
                      na.rm = TRUE)
    cor <- dat[dat$estimator == "ML", grepl("correlation", colnames(dat))]
    cor <- split(cor, rep(unique(dat$replication), each = dat$p[1] + 1))
    cor <- Reduce("+", cor) / length(cor)
    rownames(cor) <- NULL
    mcse.mse <- sqrt(colSums((abw^2 - matrix(mse, byrow = TRUE, ncol = length(mse),
                                       nrow = reps))^2, na.rm = TRUE) / (reps * (reps - 1)))
    mcse.power <- sqrt(power * (1 - power) / reps)
    return(data.frame(tmp, angle = angle, epv = epv, bias = bias, sd = sd,
                      median = median, min = min, max = max, a.bias = a.bias, 
                      rel.bias = rel.bias, mse = mse, mcse.mse = mcse.mse,
                      power = power, mcse.power = mcse.power,
                      sep = sep, converged = converged[tmp$estimator], 
                      cor))
  })
  stats <- do.call(rbind, stats)
  return(stats)
}
