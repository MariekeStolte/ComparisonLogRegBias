source("doSim.R")
sc <- data.frame(p = rep(c(1, 4, 8), times = c(4, 3, 3)), 
                 type = c("nv", "binaer.bal", "binaer.unbal", "logNv", 
                          rep(c("stetig", "binaer", "gemischt"), 2)),
                 p.nv = c(1, rep(0, 3), 2, 0, 1, 4, 0, 2),
                 p.bin.bal = c(0, 1, rep(0, 3), 2, 1, 0, 4, 2),
                 p.bin.unbal = c(0, 0, 1, 0, 0, 2, 1, 0, 4, 2),
                 p.logNv = c(rep(0, 3), 1, 2, 0, 1, 4, 0, 2))
system.time(sim.8 <- doSim(ns = c(20, 30, 50, 100, 200, 500, 1000), betas = 0:2,
                           cors = c(0, 0.8), p = 8,
                           scenarios = sc, prob.bin = 0.2,
                           mu.logNv = (-log(exp(1) - 1) - 1) / 2, sd.logNv = 1,
                           seed = 19102021, 
                           repl1 = 10, repl2 = 50))
# save(sim.8, file = "Results8.RData")
