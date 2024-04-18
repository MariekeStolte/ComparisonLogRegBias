# singleBoxplotBias1 - Creates boxplots of the bias (each n individually, without Ridge)
# Input: res.sim - (data.frame) Result from doSim
#        n, cor, scenario, beta - Parameter settings for which to plot
#                                   should
#        ylim - s.boxplot(), default: c(-5.5, 5.5)
# Output: None, boxplots split into beta and cor
#
singleBoxplotBias1 <- function(res.sim, n, cor, scenario, beta, ylim = c(-5, 5)) {
  margin <- par()$mar
  par(mar = margin + c(1, 0, -1, 0), mgp = c(3, 0.8, 0))
  on.exit(par(mar = margin))
  cols <- hcl.colors(10, "Spectral")[-(4:5)][-6]
  res.sim$estimator <- as.character(res.sim$estimator)
  res.sim <- res.sim[res.sim$estimator != "Ridge", ]
  res.sim$estimator[res.sim$estimator == "Kenne Pagui et al."] <- "KPSS"
  res.sim$estimator[res.sim$estimator == "OA"] <- "ÖA"
  res.sim$estimator <- factor(res.sim$estimator, 
                              levels = c("ML", "CM", "Haldane", "Firth", "KPSS", 
                                         "WA", "ÖA", "Separation"))
  dat <- res.sim[res.sim$n == n & res.sim$cor == cor & 
                   res.sim$scenario == scenario & abs(res.sim$beta) == beta, ]
  
  if(nrow(dat) == 0) return(invisible(NULL))
  sep <- mean(dat$separation, na.rm = TRUE)
  dat$converged[dat$estimator == "WA"] <- dat$converged[dat$estimator == "ML"]
  tmp <- dat[dat$replication == 1 & dat$estimator != "Separation", 
             c("beta", "estimator", "n", "p", "cor",  "scenario", "distr")]
  dat <- dat[!dat$separation & !is.na(dat$separation) & 
               dat$estimator != "Separation", ]
  not.bin <- !grepl("binaer", dat$scenario[1])
  if(not.bin) {
    dat <- dat[dat$estimator != "Haldane", ]
    cols <- cols[-3]
  }
  dat$estimate[!(dat$converged | is.na(dat$converged)) | abs(dat$estimate) > 10] <- NA
  dat$bias <- dat$estimate - dat$beta
  dat <- droplevels(dat)
  n <- prod(dim(with(dat, table(distr, estimator))))
  n2 <- nlevels(dat$distr)
  VEC <- (1:(n + n / n2 - 1))[-seq(n2 + 1, (n + n / n2 - 1), by = n2 + 1)]
  bp <- boxplot(bias ~ distr + estimator + n, data = dat,
                plot = FALSE, sep = ":")
  bp[["names"]] <- sapply(strsplit(bp[["names"]], ":", fixed = TRUE), 
                          function(x) x[1])
  bxp(bp, main = "", medlwd = 2,
      boxfill = rep(cols, each = length(unique(dat$distr))), 
      las = 2, cex.axis = 0.75, outcex = 0.5, xlab = "", ylab = "",
      outcol = rep(cols, each = length(unique(dat$distr))), ylim = ylim,  
      at = VEC)
  title(ylab = "Bias", line = 2, cex.lab = 0.8)
  legend("topleft", 
         legend = levels(dat$estimator),
         fill = cols, 
         xpd = TRUE, horiz = TRUE, bty = "n", inset = c(0, -0.1),
         cex = 0.75,  x.intersp = 0.5)
  legend("topright", "Arithm. Mean", pch = 19, cex = 0.75, xpd = TRUE, 
         inset = c(0, -0.1), bty = "n")
  abline(h = 0, col = "red")
  m <- aggregate(dat$bias, list(dat$distr, dat$estimator), mean, 
                 na.rm = TRUE, drop = FALSE)$x
  points(VEC, m, pch = 20)
  # print(m)
  return(invisible(NULL))
}

# singleBoxplotBias - Creates boxplots of the bias (each n individually, with Ridge)
# Input: res.sim - (data.frame) Result from doSim
#        n, cor, scenario, beta - Parameter settings for which to plot
#                                   should
#        ylim - s.boxplot(), default: c(-5.5, 5.5)
# Output: None, boxplots split into beta and cor
#
singleBoxplotBias <- function(res.sim, n, cor, scenario, beta, ylim = c(-5, 5)) {
  margin <- par()$mar
  par(mar = margin + c(1, 0, -1, 0), mgp = c(3, 0.8, 0))
  on.exit(par(mar = margin))
  cols <- hcl.colors(10, "Spectral")[-(4:5)]
  res.sim$estimator <- as.character(res.sim$estimator)
  res.sim$estimator[res.sim$estimator == "Kenne Pagui et al."] <- "KPSS"
  res.sim$estimator[res.sim$estimator == "OA"] <- "ÖA"
  res.sim$estimator <- factor(res.sim$estimator, 
                              levels = c("ML", "CM", "Haldane", "Firth", "KPSS", 
                                         "Ridge", "WA", "ÖA", "Separation"))
  dat <- res.sim[res.sim$n == n & res.sim$cor == cor & 
                   res.sim$scenario == scenario & abs(res.sim$beta) == beta, ]
  
  if(nrow(dat) == 0) return(invisible(NULL))
  # dat <- res.sim[[54]]
  sep <- mean(dat$separation, na.rm = TRUE)
  dat$converged[dat$estimator == "WA"] <- dat$converged[dat$estimator == "ML"]
  tmp <- dat[dat$replication == 1 & dat$estimator != "Separation", 
             c("beta", "estimator", "n", "p", "cor",  "scenario", "distr")]
  dat <- dat[!dat$separation & !is.na(dat$separation) & 
               dat$estimator != "Separation", ]
  not.bin <- !grepl("binaer", dat$scenario[1])
  if(not.bin) {
    dat <- dat[dat$estimator != "Haldane", ]
    cols <- cols[-3]
  }
  # set unconverged / too large values to NA
  dat$estimate[!(dat$converged | is.na(dat$converged)) | abs(dat$estimate) > 10] <- NA
  dat$bias <- dat$estimate - dat$beta
  dat <- droplevels(dat)
  n <- prod(dim(with(dat, table(distr, estimator))))
  n2 <- nlevels(dat$distr)
  VEC <- (1:(n + n / n2 - 1))[-seq(n2 + 1, (n + n / n2 - 1), by = n2 + 1)]
  bp <- boxplot(bias ~ distr + estimator + n, data = dat,
                plot = FALSE, sep = ":")
  bp[["names"]] <- sapply(strsplit(bp[["names"]], ":", fixed = TRUE), 
                          function(x) x[1])
  bxp(bp, main = "", medlwd = 2,
      boxfill = rep(cols, each = length(unique(dat$distr))),   
      las = 2, cex.axis = 0.75, outcex = 0.5, xlab = "", ylab = "",
      outcol = rep(cols, each = length(unique(dat$distr))), ylim = ylim, 
      at = VEC)
  title(ylab = "Bias", line = 2, cex.lab = 0.8)
  legend("topleft", 
         legend = levels(dat$estimator),
         fill = cols, 
         xpd = TRUE, horiz = TRUE, bty = "n", inset = c(0, -0.1),
         cex = 0.75,  x.intersp = 0.5)
  legend("topright", "Arithm. Mean", pch = 19, cex = 0.75, xpd = TRUE, 
         inset = c(0, -0.1), bty = "n")
  abline(h = 0, col = "red")
  m <- aggregate(dat$bias, list(dat$distr, dat$estimator), mean, 
                 na.rm = TRUE, drop = FALSE)$x
  points(VEC, m, pch = 20)
  return(invisible(NULL))
}

