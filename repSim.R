source("oneSim.R")

# repSim - Generates data and applies oneSim to it
# Input: x.fun - (function) Function to generate design matrix
#        beta - (numerical) true beta (appropriate dimension)
#        rep - (numeric) number of repetitions (default: 1000)
#        ... - Arguments for x.fun
# Output: data.frame with columns see oneSim() and additionally
#           n, p, cor, scenario, distr, replication
#
repSim <- function(beta, n, cor, scenarios, ind.scen, repl1 = 1000, repl2 = 5000,
                   prob.bin, mu.logNv, sd.logNv) {
  p.nv <- scenarios$p.nv[ind.scen]
  p.bin.bal <- scenarios$p.bin.bal[ind.scen]
  p.bin.unbal <- scenarios$p.bin.unbal[ind.scen]
  p.logNv <- scenarios$p.logNv[ind.scen]
  problem <- grepl("binaer", scenarios$type[ind.scen]) ||
    (scenarios$p == 8 && cor != 0 && beta != 0 && scenarios$type == "gemischt") ||
    (n == 20 && scenarios$p == 8 && cor == 0 && beta == 0 && scenarios$type == "gemischt") ||
    (n == 20 && scenarios$p == 4 && cor == 0 && beta == 2 && scenarios$type == "gemischt") ||
    (n < 200 && scenarios$p == 8 && cor == 0 && beta > 0) ||
    (n < 200 && scenarios$p > 1 && cor == 0.8 && beta > 0)
  reps <- ifelse(problem, repl2, repl1)
  res <- replicate(reps, {
    tmp <- oneSim(n = n, p.nv = p.nv, p.bin.bal = p.bin.bal, 
                  p.bin.unbal = p.bin.unbal, p.logNv = p.logNv, prob.bin = prob.bin,
                  cor = cor, beta = beta, mu.logNv = mu.logNv, sd.logNv = sd.logNv)
  }, simplify = FALSE) 
  res <- do.call(rbind, res)
  res <- data.frame(res, n = n, p = scenarios$p[ind.scen], cor = cor, 
                    scenario = scenarios$type[ind.scen],
                    distr = factor(c("Intercept", rep("N(0,1)", p.nv),  
                                     rep("Ber(0.5)", p.bin.bal),
                                     rep(paste0("Ber(", prob.bin, ")"), p.bin.unbal),
                                     rep(paste0("logN(", round(mu.logNv, 2), ", ",
                                                round(sd.logNv, 2), ")"), p.logNv))),
                    replication = rep(1:reps, each = nrow(res) / reps))

  # save(res, file = paste0("res_", ind.scen, "_", 
  #                         gsub(":", "_", Sys.time(), fixed = TRUE), ".RData"))
  return(res)
}

