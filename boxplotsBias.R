# boxplotsBias - Creates boxplots of the bias (all ns in one plot)
# Input: res.sim - (data.frame) Result from doSim
#        ylim - s.boxplot(), default: c(-5.5, 5.5)
# Output: None, boxplots split into beta and cor
#
boxplotsBias <- function(res.sim, ylim = c(-5.5, 5.5)) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mar = op$mar - c(1, 0, 2, 0), mgp = c(3, 0.9, 0))
  cols <- hcl.colors(10, "Spectral")[-(4:5)][-6]
  res.sim$estimator <- as.character(res.sim$estimator)
  # levels(res.sim$estimator) <- c(levels(res.sim$estimator), "KPSS", "\\\"OA")
  res.sim$estimator[res.sim$estimator == "Kenne Pagui et al."] <- "KPSS"
  res.sim$estimator[res.sim$estimator == "OA"] <- "ÖA"
  res.sim$estimator <- factor(res.sim$estimator, 
                              levels = c("ML", "CM", "Haldane", "Firth", "KPSS", 
                                         "WA", "ÖA", "Separation"))
  res.sim <- split(res.sim, list(res.sim$scenario, res.sim$cor, abs(res.sim$beta)), 
                   drop = TRUE)
  stats <- lapply(res.sim, function(dat) {
    if(nrow(dat) == 0) return(invisible(NULL))
    sep <- aggregate(dat$separation, list(dat$n), mean, 
                     na.rm = TRUE, drop = FALSE)$x
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
    n <- prod(dim(with(dat, table(distr, estimator, n))))
    n2 <- prod(dim(with(dat, table(distr, estimator))))
    VEC <- (1:(n + n / n2 - 1))[-seq(n2 + 1, n, by = n2 + 1)]
    bp <- boxplot(bias ~ distr + estimator + n, data = dat,
                  plot = FALSE, sep = ":")
    names <- sapply(strsplit(bp[["names"]], ":", fixed = TRUE), 
                    function(x) x[1])
    bp[["names"]] <- ""
    bxp(bp, 
        medlwd = 2,
        boxfill = rep(cols, each = length(unique(dat$distr))), 
        xlim = c(ifelse(any(sep == 1), 0, length(unique(dat$distr)) + 1), 
                 ifelse(any(sep == 1), max(VEC), max(VEC) - length(unique(dat$distr)))),  
        las = 2, cex.axis = 0.8, outcex = 0.3, xlab = "", ylab = "",
        outcol = rep(cols, each = length(unique(dat$distr))), ylim = ylim, 
        show.names = FALSE, at = VEC)
    ins <- switch(as.character(max(VEC)), "104" = 5.5, "90" = 4.75, "153" = 7.65, 
                  "216" = 10.75, "132" = 6.55, "131" = 6.9, "185" = 9.7, 
                  "43" = 2.25, "154" = 8.1)
    text(par("usr")[1] - ins, 0, "Bias",
         cex = 0.8, srt = 270, xpd = TRUE)
    legend("topleft", 
           legend = levels(dat$estimator), 
           fill = cols,
           xpd = TRUE, horiz = TRUE, bty = "n", inset = c(0, -0.05),
           cex = 0.65,  x.intersp = 0.5)
    legend("topright", "Mean", pch = 19, cex = 0.75, xpd = TRUE, 
           inset = c(0, -0.05), bty = "n")
    pos <- seq(length(VEC) - 3 * length(unique(dat$distr)) - 1, length(VEC), by = 3)
    text(max(VEC), ylim[1],
         paste0(c("Coefficient:", names[length(unique(dat$distr)):1]), collapse = "\n"),
         cex = 0.6, adj = c(1, 1), font = 3, srt = 270)
    mtext(side = 1, line = 0.5, at = colMeans(matrix(VEC, ncol = n / n2)),
          text = do.call(expression, sapply(unique(dat$n), function(m) bquote(italic(n) == .(m)))), 
          cex = 0.75)
    mtext(side = 1, line = 1.5, at = colMeans(matrix(VEC, ncol = n / n2)),
          text = paste0("(Separable: ",
                        format(sep[sep != 1], digits = 3), ")"),
          cex = 0.75)
    m <- aggregate(dat$bias, list(dat$distr, dat$estimator, dat$n), mean, 
                   na.rm = TRUE, drop = FALSE)
    m <- m$x
    points(VEC, m, pch = 20, cex = 0.6)
  })
  return(invisible(NULL))
}


# boxplotsBias2 - Creates boxplots of the bias (all ns in one plot) incl. Ridge estimator
# Input: res.sim - (data.frame) Result from doSim
#        ylim - s.boxplot(), default: c(-5.5, 5.5)
# Output: None, boxplots split into beta and cor
#

boxplotsBias2 <- function(res.sim, ylim = c(-5.5, 5.5)) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mar = op$mar - c(1, 0, 2, 0), mgp = c(3, 0.9, 0))
  cols <- hcl.colors(10, "Spectral")[-(4:5)]
  res.sim$estimator <- as.character(res.sim$estimator)
  res.sim$estimator[res.sim$estimator == "Kenne Pagui et al."] <- "KPSS"
  res.sim$estimator[res.sim$estimator == "OA"] <- "ÖA"
  res.sim$estimator <- factor(res.sim$estimator, 
                              levels = c("ML", "CM", "Haldane", "Firth", "KPSS", 
                                         "Ridge", "WA", "ÖA", "Separation"))
  res.sim <- split(res.sim, list(res.sim$scenario, res.sim$cor, abs(res.sim$beta)), 
                   drop = TRUE)
  stats <- lapply(res.sim, function(dat) {
    if(nrow(dat) == 0) return(invisible(NULL))
    sep <- aggregate(dat$separation, list(dat$n), mean, 
                     na.rm = TRUE, drop = FALSE)$x
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
    n <- prod(dim(with(dat, table(distr, estimator, n))))
    n2 <- prod(dim(with(dat, table(distr, estimator))))
    if(n2 == n) {
      VEC <- (1:(n + n / n2 - 1))
    } else {
      VEC <- (1:(n + n / n2 - 1))[-seq(n2 + 1, n, by = n2 + 1)]
    }
    bp <- boxplot(bias ~ distr + estimator + n, data = dat,
                  plot = FALSE, sep = ":")
    names <- sapply(strsplit(bp[["names"]], ":", fixed = TRUE), 
                    function(x) x[1])
    bp[["names"]] <- ""
    bxp(bp, 
        medlwd = 2,
        boxfill = rep(cols, each = length(unique(dat$distr))), 
        xlim = c(ifelse(any(sep == 1), 0, length(unique(dat$distr)) + 1), 
                 ifelse(any(sep == 1), max(VEC), max(VEC) - length(unique(dat$distr)))),  
        las = 2, cex.axis = 0.8, outcex = 0.3, xlab = "", ylab = "",
        outcol = rep(cols, each = length(unique(dat$distr))), ylim = ylim, 
        show.names = FALSE, at = VEC)
    ins <- switch(as.character(max(VEC)), "104" = 5.5, "90" = 4.75, "153" = 7.65, 
                  "216" = 10.75, "132" = 6.55, "131" = 6.9, "185" = 9.7, 
                  "43" = 2.25, "154" = 8.1)
    text(par("usr")[1] - ins, 0, "Bias",
         cex = 0.8, srt = 270, xpd = TRUE)
    legend("topleft", 
           legend = levels(dat$estimator), 
           fill = cols,
           xpd = TRUE, horiz = TRUE, bty = "n", inset = c(0, -0.05),
           cex = 0.65,  x.intersp = 0.5)
    legend("topright", "Mean", pch = 19, cex = 0.75, xpd = TRUE, 
           inset = c(0, -0.05), bty = "n")
    pos <- seq(length(VEC) - 3 * length(unique(dat$distr)) - 1, length(VEC), by = 3)
    text(max(VEC), ylim[1],
         paste0(c("Coefficient:", names[length(unique(dat$distr)):1]), collapse = "\n"),
         cex = 0.6, adj = c(1, 1), font = 3, srt = 270)
    abline(h = 0, col = "red")
    mtext(side = 1, line = 0.5, at = colMeans(matrix(VEC, ncol = n / n2)),
          text = do.call(expression, sapply(unique(dat$n), function(m) bquote(italic(n) == .(m)))), 
          cex = 0.75)
    mtext(side = 1, line = 1.5, at = colMeans(matrix(VEC, ncol = n / n2)),
          text = paste0("(Separable: ",
                        format(sep[sep != 1], digits = 3), ")"),
          cex = 0.75)
    m <- aggregate(dat$bias, list(dat$distr, dat$estimator, dat$n), mean, 
                   na.rm = TRUE, drop = FALSE)
    m <- m$x
    points(VEC, m, pch = 20, cex = 0.6)
  })
  return(invisible(NULL))
}

load("Results1.RData", verbose = TRUE)
sim.1 <- sim.1[sim.1$estimator != "Ridge", ] # all NAs
# pdf("../Plots/BoxplotBias1_en.pdf", width = 9.7, height = 6.2)
boxplotsBias(sim.1)
# dev.off()



load("Results4.RData", verbose = TRUE)
# pdf("../Plots/BoxplotBias4_en.pdf", width = 9.7, height = 6.2)
boxplotsBias2(sim.4)
# dev.off()
rm(sim.1, sim.4)


load("Results8.RData", verbose = TRUE)
# pdf("../Plots/BoxplotBias8_en.pdf", width = 9.7, height = 6.2)
boxplotsBias2(sim.8)
# dev.off()
rm(sim.8)
