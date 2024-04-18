# plotStats - Plots results
# Input: stats - (data.frame) Result from calcStats
#        what - (character) Column name, what should be plotted?
#        p - (numeric) For which p should results be plotted?
#        ylab - (character) ylab for plots (default: what)
#        cols - (character or numeric) color vector
#           (default: hcl.colors(9, "Spectral")[-(4:5)][c(1:2, 7, 6, 4, 3, 5)])
#        est.names - (character) Names of the estimators for whom the results are to be returned
#           should be plotted (default: unique(droplevels(stats$estimator)))
#        names.legend - (named character) guess names for the legend
#           (default: c(ML = "ML", CM = "CM", OA = "ÖA",
#           WA = "WA", Firth = "Firth", Haldane = "Haldane",
#           “Kenne Pagui et al.” = "KPSS")), names(names.legend) must be the 
#           names from est.names
#        max.ylim - (numeric) 2nd entry for ylim (default: 1)
#        min.ylim - (numeric) 1st entry for ylim (default: 0)
#        inset.y - (numeric) y-inset for legend (default: -0.2)
#        inset.x - (numeric) x-inset for legend (default: rep(0, 5)), i-th
#             Entry is inset at i columns in the plot
# Output: Plotted menus as a table
#
plotStats <- function(stats, what, p, ylab = what, 
                      cols = hcl.colors(9, "Spectral")[-(4:5)][c(1:2, 7, 6, 4, 3, 5)],
                      est.names = unique(droplevels(stats$estimator)), 
                      names.legend = c(ML = "ML", CM = "CM", OA = "ÖA", 
                                       WA = "WA", Firth = "Firth", Haldane = "Haldane",
                                       "Kenne Pagui et al." = "KPSS", Ridge = "Ridge"),
                      max.ylim = 1, 
                      min.ylim = 0, inset.y = -0.2, inset.x = rep(0, 5)) {
  stats <- droplevels(stats)
  stats <- stats[stats$p == p, c("beta", "estimator", "n", "p", "cor", "scenario",
                                 "distr", "angle", "epv", "sep", what, paste0("mcse.", what))]
  xlim <- c(18, 1200)
  names(cols) <- est.names
  s <- split(stats, f = list(stats$distr, abs(stats$beta), stats$cor, stats$scenario), 
             drop = TRUE)
  means <- lapply(s, function(m) {
    aggregate(m[, what], list(m$estimator, m$n), mean, na.rm = TRUE)
  })
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  means <- lapply(means, function(m) {
    colnames(m) <- c("estimator", "n", what)
    return(m)
  })
  ses <- lapply(s, function(m) {
    aggregate(m[, paste0("mcse.", what)], list(m$estimator, m$n), function(x) x)
  })
  ses <- lapply(ses, function(m) {
    colnames(m) <- c("estimator", "n", paste0("mcse.", what))
    return(m)
  })
  if(p == 1) {
    par(mar = c(0,ifelse(what == "power", 1.5, 1), 0.5,0), oma=c(7, 4, 2.5, 0.5), 
        mgp = c(1, 0.9, 0))
    plot.rows <- ifelse(what == "rel.bias", 2, 3)
    plot.cols <- 2
    par(mfrow = c(plot.rows, 2))
    counter <- 0
    if(what == "power") {
      ylims <- rep(1, length(means))
    } else {
      ylims <- sapply(1:length(means), function(i) max(na.omit(as.numeric(means[[i]][, what] + 
                                                                            ses[[i]][, paste0("mcse.", what)]))))
      ylims <- ifelse(is.finite(ylims), ylims, max.ylim)
    }
    ylim <- max(ylims[1:(plot.rows * 2)])
    for(i in seq(along = means)) {
      if(counter == 2 * plot.rows) {
        legend("bottomleft", names.legend[as.character(est.names)], col = cols, lty = 1, pch = 16, 
               inset = c(inset.x[plot.cols], inset.y), # cex = 0.9,
               xpd = NA, bty = "n", title = "Estimator:", horiz = TRUE)
        ylim <- max(ylims[i:(plot.rows * 2 + i - 1)])
        counter <- 0
      }
      if(what == "rel.bias" && s[[i]]$beta[1] == 0) next
      m <- means[[i]]
      se <- ses[[i]]
      y.lab <- ifelse(what == "power" & s[[i]]$beta[1] == 0, "Type I Error Rate", ylab)
      plot(NA, xlim = xlim, xlab = "n", log = "x",
           ylim = c(min.ylim, ifelse(what == "power" & s[[i]]$beta[1] == 0, 0.1, ylim)),
           ylab = "", xaxt = "n",
           yaxt = ifelse(counter %% 2 == 0, "s", "n"), main = "", las = 1)
      axis(1, at = s[[i]]$n, 
           labels = ifelse(rep(counter >= (plot.rows - 1) * 2, length(s[[i]]$n)), s[[i]]$n, ""), 
           las = 1)
      if(counter %% 2 == 0) {
        mtext(side = 2, text = y.lab, line = ifelse(what == "power", 2.8, 2.2), cex = 0.8)
        mtext(side = 2, text = bquote(beta == .(abs(s[[i]]$beta[1]))), 
              line = ifelse(what == "power", 4.1, 3.5))
      }
      if(counter >= (plot.rows - 1) * 2) {
        mtext(side = 1, text = "n", line = 2.8, cex = 0.9)
      }
      mtext(side = 3, text = ifelse(counter < 2, as.character(s[[i]]$distr[1]), ""), 
            line = 0.8, cex = 0.9)
      abline(h = 0, lty = 2)
      if(what == "power") abline(h = 1, lty = 2)
      if(what == "power" & abs(s[[i]]$beta[1]) == 0) abline(h = 0.05, col = "red", lty = 2)
      for(est in est.names) {
        l <- m[m$estimator == est, ]
        lines(l$n, l[, what], col = cols[est], type = "b", pch = 16)
        s.e <- se[se$estimator == est, ]
        arrows(x0 = s.e$n, y0 = l[, what], 
               y1 = l[, what] + s.e[, paste0("mcse.", what)], angle = 90, 
               col = cols[est], length = 0.05)
        arrows(x0 = s.e$n, y0 = l[, what] , 
               y1 = l[, what] - s.e[, paste0("mcse.", what)], angle = 90, 
               col = cols[est], length = 0.05)
      }
      counter <- counter + 1
    }
    leg <- as.character(est.names)
    ins <- inset.x[plot.cols]
    cl <- cols
    if(what == "mse" && s[[i]]$scenario[1] != "binaer") {
      leg <- leg[leg != "Haldane"]
      ins <- ins + ifelse(s[[i]]$scenario[1] == "stetig", 0.33, 0.55)
      cl <- cl[names(cl) != "Haldane"]
    }
    if(what == "power") {
      leg <- leg[!leg %in% c("Haldane", "OA", "WA", "Ridge")]
      cl <- cl[!names(cl) %in% c("Haldane", "OA", "WA", "Ridge")]
    }
    legend("bottomleft", names.legend[leg], col = cl, lty = 1, pch = 16, 
           inset = c(ins, inset.y), 
           xpd = NA, bty = "n", title = "Estimator:", horiz = TRUE)
    
  } else {
    par(mar=c(0,ifelse(what == "power",1.5, 1),0.5,0), oma=c(7,4,2.5,0.5), mgp = c(1, 0.9, 0))
    plot.cols <- length(unique(stats$distr[stats$scenario == s[[1]]$scenario[1]]))
    plot.rows <- ifelse(what == "rel.bias", 2, 3) # 1
    counter <- 0
    par(mfrow = c(plot.rows, plot.cols))
    if(what == "power") {
      ylims <- rep(1, length(means))
    } else {
      ylims <- sapply(means, function(m) max(na.omit(as.numeric(m[, what]))))
      ylims <- ifelse(is.finite(ylims), ylims, max.ylim)
    }
    ylim <- max(ylims[1:(plot.rows * plot.cols)])
    for(i in seq(along = means)) {
      if(counter == plot.cols * plot.rows) {
        leg <- as.character(est.names)
        ins <- inset.x[plot.cols]
        cl <- cols
        if(what == "mse" && s[[i-1]]$scenario[1] != "binaer") {
          leg <- leg[leg != "Haldane"]
          cl <- cl[names(cl) != "Haldane"]
          ins <- ins + ifelse(s[[i-1]]$scenario[1] == "stetig", 0.33, 0.55)
        }
        if(what == "power") {
          leg <- leg[!leg %in% c("Haldane", "OA", "WA", "Ridge")]
          cl <- cl[!names(cl) %in% c("Haldane", "OA", "WA", "Ridge")]
        }
        legend("bottomleft", names.legend[leg], col = cl, lty = 1, pch = 16, 
               inset = c(ins, inset.y),
               xpd = NA, bty = "n", title = "Estimator:", horiz = TRUE)
        counter <- 0
        plot.cols <- length(unique(stats$distr[stats$scenario == s[[i]]$scenario[1]]))
        par(mfrow = c(plot.rows, plot.cols))
        ylim <- max(ylims[i:(plot.rows * plot.cols + i - 1)])
      }
      if(what == "rel.bias" && s[[i]]$beta[1] == 0) next
      m <- means[[i]]
      se <- ses[[i]]
      y.lab <- ifelse(what == "power" & s[[i]]$beta[1] == 0, "Type I Error Rate", ylab)
      plot(NA, xlim = xlim, xlab = "$n$", log = "x",
           ylim = c(min.ylim, ifelse(what == "power" & s[[i]]$beta[1] == 0, 0.1, ylim)),
           ylab = "", xaxt = "n",
           yaxt = ifelse(counter %% plot.cols == 0, "s", "n"), main = "", las = 1)
      axis(1, at = s[[i]]$n, 
           labels = ifelse(rep(counter >= (plot.rows - 1) * plot.cols, length(s[[i]]$n)), s[[i]]$n, ""), las = 3)
      if(counter %% plot.cols == 0) {
        mtext(side = 2, text = y.lab, line = ifelse(what == "power", 2.8, 2.2), cex = 0.8)
        mtext(side = 2, text = bquote(beta == .(abs(s[[i]]$beta[1]))), line = ifelse(what == "power", 4.1, 3.5))
      }
      if(counter >= (plot.rows - 1) * plot.cols) {
        mtext(side = 1, text = "n", line = 2.8, cex = 0.9)
      }
      mtext(side = 3, text = ifelse(counter < plot.cols, as.character(s[[i]]$distr[1]), ""), 
            line = 0.8, cex = 0.9)
      abline(h = 0, lty = 2)
      if(what == "power") abline(h = 1, lty = 2)
      if(what == "power" & abs(s[[i]]$beta[1]) == 0) abline(h = 0.05, col = "red", lty = 2)
      for(est in est.names) {
        l <- m[m$estimator == est, ]
        lines(l$n, l[, what], col = cols[est], type = "b", pch = 16)
        s.e <- se[se$estimator == est, ]
        arrows(x0 = s.e$n, y0 = l[, what], 
               y1 = l[, what] + s.e[, paste0("mcse.", what)], angle = 90, 
               col = cols[est], length = 0.05)
        arrows(x0 = s.e$n, y0 = l[, what] , 
               y1 = l[, what] - s.e[, paste0("mcse.", what)], angle = 90, 
               col = cols[est], length = 0.05)
      }
      counter <- counter + 1
    }
    leg <- as.character(est.names)
    ins <- inset.x[plot.cols]
    cl <- cols
    if(what == "mse" && s[[i]]$scenario[1] != "binaer") {
      leg <- leg[leg != "Haldane"]
      ins <- ins + ifelse(s[[i]]$scenario[1] == "stetig", 0.33, 0.55)
      cl <- cl[names(cl) != "Haldane"]
    }
    if(what == "power") {
      leg <- leg[!leg %in% c("Haldane", "OA", "WA", "Ridge")]
      cl <- cl[!names(cl) %in% c("Haldane", "OA", "WA", "Ridge")]
    }
    legend("bottomleft", names.legend[leg], col = cl, lty = 1, pch = 16, 
           inset = c(ins, inset.y), # cex = 0.9,
           xpd = NA, bty = "n", title = "Estimator:", horiz = TRUE)
  }
  return(invisible(means))
}

getCaptionsStats <- function(stats, what) {
  stats <- droplevels(stats)
  s <- split(stats, f = list(stats$cor, stats$scenario), 
             drop = TRUE)
  cap <- numeric(length(s))
  for(i in seq(along = s)){
    if(s[[i]]$p[1] == 1) {
      s[[i]]$dist <- sapply(s[[i]]$scenario, function(s)
        switch(as.character(s), "binaer.bal" = "$Ber(0.5)$",
               "binaer.unbal" = "$Ber(0.2)$",
               "logNv" = "$logN(-0.77, 1)$",
               "nv" = "$N(0,1)$"))
    } else {
      s[[i]]$dist <- gsub("binaer", "binary", s[[i]]$scenario)
      s[[i]]$dist <- gsub("stetig", "continuous", s[[i]]$scenario)
      s[[i]]$dist <- gsub("gemischt", "mixed", s[[i]]$scenario)
    }
    cap[i] <- paste0(what, " for ", s[[i]]$dist[1],
                      ifelse(s[[i]]$p[1] == 1, " distributed ", ""), "data, ",
                      "$p = ", s[[i]]$p[1], "$ und $\\rho = ", s[[i]]$cor[1], 
                      "$")
  }
  return(cap)
}
