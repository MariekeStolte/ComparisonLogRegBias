# singlePlotDistr - Bias separated by total, overlap and separable
# (for sim.4, n = 20, beta = 1, mixed, rho = 0)
# Input: est - (character) Name of the estimator for which the value should be evaluated
# Output: Plot
#
singlePlotDistr <- function(est){
  op <- par(no.readonly = TRUE)
  par(mgp = c(1.5, 0.9, 0), mar = op$mar + c(0, 0, 0, 2.5))
  ylim <- c(-5.5, 5.5)
  res.sim <- sim.4[sim.4$n == 20 & abs(sim.4$beta) == 1 & sim.4$scenario == "gemischt" & sim.4$cor == 0, ]
  cols <- hcl.colors(6, "Spectral")[6:4]
  res.sim <- res.sim[res.sim$estimator == est,]
  res.sim <- droplevels(res.sim)
  res.sim$estimate[(!(res.sim$converged | is.na(res.sim$converged))) | abs(res.sim$estimate) > 10] <- NA
  res.sim$bias <- res.sim$estimate - res.sim$beta
  s <- split(res.sim, f = list(abs(res.sim$beta), res.sim$cor, res.sim$scenario))
  i <- 1
  a <- rep(s[[i]]$bias, 2)
  sep <- c(c("Overlap", "Separable")[s[[i]]$separation + 1], 
           rep("Overall", nrow(s[[i]])))
  sep <- factor(sep, levels =  c("Overall", "Overlap", "Separable"))
  s[[i]] <- droplevels(s[[i]])
  bp <- boxplot(a ~ rep(s[[i]]$distr, 2) + sep +  rep(s[[i]]$n, 2), plot = FALSE, sep = ":")
  bp[["names"]] <- sapply(strsplit(bp[["names"]], ":", fixed = TRUE), 
                          function(x) x[1])
  n <- 15
  n2 <- 5
  VEC <- (1:(n + n / n2 - 1))[-seq(n2 + 1, n, by = n2 + 1)]
  bxp(bp, main = "",
      medlwd = 2,
      boxfill = rep(cols, each = length(unique(s[[i]]$distr))), #xlim = c(2, max(VEC) - 1),  
      las = 2, cex.axis = 0.75, outcex = 0.5, xlab = "", ylab = "Bias",
      outcol = rep(cols, each = length(unique(s[[i]]$distr))), ylim = ylim, at = VEC) 
  legend("topright", legend = c(sort(unique(as.character(sep))), "Arithm. Mean"), 
         fill = c(cols, "white"), border = c(rep("black", 3), "white"), 
         bty = "n", pch = c(rep(NA, 3), 20),
         cex = 0.60, xpd = TRUE, inset = c(-0.2, 0))
  abline(h = 0, col = "red")
  s[[i]]$separation <- factor(s[[i]]$separation, levels = c(FALSE, TRUE))
  prop <- format(as.numeric(prop.table(table(s[[i]]$separation, s[[i]]$n), 2)), digits = 4)
  mtext("Proportion:", at = 3, cex = 0.75)
  mtext(prop, at = c(9, 15), cex = 0.7)
  m <- aggregate(a, list(rep(s[[i]]$distr, 2), sep, rep(s[[i]]$n, 2)), mean, 
                 na.rm = TRUE, drop = FALSE)$x
  points(VEC, m, cex = 0.6, pch = 20)
  par(op)
}

load("Results4.RData")
# pdf("../Plots/plotBiasSep.pdf", height = 4, width = 6)
singlePlotDistr("Firth")
singlePlotDistr("Kenne Pagui et al.")
# dev.off()