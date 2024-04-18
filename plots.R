source("plotStats.R")
source("singleBoxplotBias.R")
source("plotH0.R")
load("Results1.RData", verbose = TRUE)
load("Results4.RData", verbose = TRUE)
load("Stats1.RData", verbose = TRUE)
load("Stats4.RData", verbose = TRUE)
load("Stats8.RData", verbose = TRUE)
################################################################################
#                               Chosen Plots                                   #
#                                                                              #
################################################################################
# Bias:
# pdf("../Plots/boxplotBiasPaper.pdf", height = 5, width = 7)
singleBoxplotBias(sim.4, n = 20, beta = 0, cor = 0, scenario = "gemischt")
singleBoxplotBias(sim.4, n = 20, beta = 1, cor = 0, scenario = "gemischt")
singleBoxplotBias1(sim.1, n = 20, beta = 1, cor = 0, scenario = "binaer.bal")
singleBoxplotBias1(sim.1, n = 20, beta = 1, cor = 0, scenario = "logNv")
singleBoxplotBias(sim.4, n = 20, beta = 2, cor = 0, scenario = "gemischt")
singleBoxplotBias(sim.4, n = 100, beta = 2, cor = 0, scenario = "gemischt")
# dev.off()


# MSE:
# pdf("../Plots/plotMSEPaper.pdf", height = 6, width = 7)
stats1 <- stats1[stats1$estimator != "Ridge", ]
plotStats(stats1, what = "mse", p = 1, max.ylim = 2, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)][-6], 
          inset.x = c(0, -0.85), inset.y = -0.6)
stats.red <- stats4[stats4$scenario == "gemischt" & stats4$cor == 0, ]
plotStats(stats.red, what = "mse", p = 4, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-4.2, 5), inset.y = -0.65)

stats.red <- stats4[stats4$scenario == "binaer" & stats4$cor == 0, ]
plotStats(stats.red, what = "mse", p = 4, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-2.05, 5), inset.y = -0.65)

stats.red <- stats4[stats4$scenario == "stetig" & stats4$cor == 0, ]
plotStats(stats.red, what = "mse", p = 4, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-2, 5), inset.y = -0.65)


stats.red <- stats8[stats8$scenario == "gemischt" & stats8$cor == 0, ]
plotStats(stats.red, what = "mse", p = 8, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-4.2, 5), inset.y = -0.65)#, 
stats.red <- stats8[stats8$scenario == "binaer" & stats8$cor == 0, ]
plotStats(stats.red, what = "mse", p = 8, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-2.05, 5), inset.y = -0.65)

stats.red <- stats8[stats8$scenario == "stetig" & stats8$cor == 0, ]
plotStats(stats.red, what = "mse", p = 8, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-2, 5), inset.y = -0.65)

stats.red <- stats8[stats8$scenario == "gemischt" & stats8$cor == 0.8, ]
plotStats(stats.red, what = "mse", p = 8, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-4.2, 5), inset.y = -0.65)#, 
stats.red <- stats8[stats8$scenario == "binaer" & stats8$cor == 0.8, ]
plotStats(stats.red, what = "mse", p = 8, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-2.05, 5), inset.y = -0.65)

stats.red <- stats8[stats8$scenario == "stetig" & stats8$cor == 0.8, ]
plotStats(stats.red, what = "mse", p = 8, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-2, 5), inset.y = -0.65)

stats.red <- stats4[stats4$scenario == "gemischt" & stats4$cor == 0.8, ]
plotStats(stats.red, what = "mse", p = 4, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-4.2, 5), inset.y = -0.65)
stats.red <- stats4[stats4$scenario == "binaer" & stats4$cor == 0.8, ]
plotStats(stats.red, what = "mse", p = 4, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-2.05, 5), inset.y = -0.65)

stats.red <- stats4[stats4$scenario == "stetig" & stats4$cor == 0.8, ]
plotStats(stats.red, what = "mse", p = 4, max.ylim = 8, ylab = "MSE", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-2, 5), inset.y = -0.65)
# dev.off()

# Tests:
# p-values under H0:
sim.red <- sim.4[sim.4$n %in% c(20, 500) & sim.4$scenario == "gemischt" & sim.4$cor == 0, ]
par(mfrow = c(2, 2))
plotH0(sim.red, ylim = c(0, 1.4))
par(mfrow = c(1, 1))
# Power:
# pdf("../Plots/plotPowerPaper.pdf", height = 6, width = 7)
plotStats(stats4[stats4$scenario == "gemischt" & stats4$cor == 0, ], what = "power", 
          p = 4, max.ylim = 1, ylab = "Power",
          names.legend = c(ML = "ML", CM = "CM", 
                           Firth = "Firth", 
                           "Kenne Pagui et al." = "KPSS"),
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-3, 5), 
          inset.y = -0.65)
plotStats(stats4[stats4$scenario == "gemischt" & stats4$cor == 0.8, ], what = "power", 
          p = 4, max.ylim = 1, ylab = "Power", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], inset.x = rep(-3, 5), 
          inset.y = -0.65)


plotStats(stats1[stats1$scenario == "binaer.bal", ], what = "power", 
          p = 1, max.ylim = 1, ylab = "Power",
          names.legend = c(ML = "ML", CM = "CM",
                           Firth = "Firth", 
                           "Kenne Pagui et al." = "KPSS"),
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-0.5, 5), 
          inset.y = -0.65)
plotStats(stats1[stats1$scenario == "binaer.unbal", ], what = "power", 
          p = 1, max.ylim = 1, ylab = "Power", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], inset.x = rep(-0.5, 5), 
          inset.y = -0.65)
plotStats(stats1[stats1$scenario == "nv", ], what = "power", 
          p = 1, max.ylim = 1, ylab = "Power",
          names.legend = c(ML = "ML", CM = "CM", 
                           Firth = "Firth", 
                           "Kenne Pagui et al." = "KPSS"),
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-0.5, 5), 
          inset.y = -0.65)
plotStats(stats1[stats1$scenario == "logNv", ], what = "power", 
          p = 1, max.ylim = 1, ylab = "Power", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], inset.x = rep(-0.5, 5), 
          inset.y = -0.65)

plotStats(stats4[stats4$scenario == "binaer" & stats4$cor == 0, ], what = "power", 
          p = 4, max.ylim = 1, ylab = "Power",
          names.legend = c(ML = "ML", CM = "CM", 
                           Firth = "Firth", 
                           "Kenne Pagui et al." = "KPSS"),
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-1.25, 5), 
          inset.y = -0.65)
plotStats(stats4[stats4$scenario == "binaer" & stats4$cor == 0.8, ], what = "power", 
          p = 4, max.ylim = 1, ylab = "Power", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], inset.x = rep(-1.25, 5), 
          inset.y = -0.65)


plotStats(stats4[stats4$scenario == "stetig" & stats4$cor == 0, ], what = "power", 
          p = 4, max.ylim = 1, ylab = "Power",
          names.legend = c(ML = "ML", CM = "CM", 
                           Firth = "Firth",
                           "Kenne Pagui et al." = "KPSS"),
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-1.25, 5), 
          inset.y = -0.65)
plotStats(stats4[stats4$scenario == "stetig" & stats4$cor == 0.8, ], what = "power", 
          p = 4, max.ylim = 1, ylab = "Power", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], inset.x = rep(-1.25, 5), 
          inset.y = -0.65)

plotStats(stats8[stats8$scenario == "gemischt" & stats8$cor == 0, ], what = "power", 
          p = 8, max.ylim = 1, ylab = "Power",
          names.legend = c(ML = "ML", CM = "CM", 
                           Firth = "Firth",
                           "Kenne Pagui et al." = "KPSS"),
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-3, 5), 
          inset.y = -0.65)
plotStats(stats8[stats8$scenario == "gemischt" & stats8$cor == 0.8, ], what = "power", 
          p = 8, max.ylim = 1, ylab = "Power", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], inset.x = rep(-3, 5), 
          inset.y = -0.65)

plotStats(stats8[stats8$scenario == "binaer" & stats8$cor == 0, ], what = "power", 
          p = 8, max.ylim = 1, ylab = "Power",
          names.legend = c(ML = "ML", CM = "CM", 
                           Firth = "Firth",
                           "Kenne Pagui et al." = "KPSS"),
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-1.25, 5), 
          inset.y = -0.65)
plotStats(stats8[stats8$scenario == "binaer" & stats8$cor == 0.8, ], what = "power", 
          p = 8, max.ylim = 1, ylab = "Power", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], inset.x = rep(-1.25, 5), 
          inset.y = -0.65)


plotStats(stats8[stats8$scenario == "stetig" & stats8$cor == 0, ], what = "power", 
          p = 8, max.ylim = 1, ylab = "Power",
          names.legend = c(ML = "ML", CM = "CM",
                           Firth = "Firth", 
                           "Kenne Pagui et al." = "KPSS"),
          cols =  hcl.colors(10, "Spectral")[-(4:5)], 
          inset.x = rep(-1.25, 5), 
          inset.y = -0.65)
plotStats(stats8[stats8$scenario == "stetig" & stats8$cor == 0.8, ], what = "power", 
          p = 8, max.ylim = 1, ylab = "Power", 
          names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                           WA = "WA", Firth = "Firth", Haldane = "Haldane",
                           "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"), 
          cols =  hcl.colors(10, "Spectral")[-(4:5)], inset.x = rep(-1.25, 5), 
          inset.y = -0.65)
# dev.off()
