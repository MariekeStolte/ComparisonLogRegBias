library(brglm2)
library(detectseparation)
library(glmnet)
source("Oezkale_Arican_2018.R")
source("generateData.R")
source("groupData.R")

# oneSim - Performs simulation on a data set
# Input: data - (data.frame or numeric matrix) data set with 1 column y
#        beta - (numeric vector) True beta
# Output: (data.frame) with columns estimate, p.value, beta, epv, angle,
#           separation (TRUE/FALSE), converged (TRUE/FALSE), estimator,
#           trial (how often did data have to be redrawn?), correlation.Xj
#
oneSim <- function(n, p.nv, p.bin.bal, p.bin.unbal, p.logNv, cor, beta, 
                   mu.logNv, sd.logNv, prob.bin) {
  for(i in 1:20) {
    data <- generateData(n = n, p.nv = p.nv, p.bin.bal = p.bin.bal, 
                         p.bin.unbal = p.bin.unbal, p.logNv = p.logNv, 
                         prob.bin = prob.bin, cor = cor, beta, mu.nv = 0, 
                         mu.logNv = mu.logNv, sd.nv = 1, sd.logNv = sd.logNv)
    if(!any(colSums(data) %in% c(0, n))) break
  }
  p <- p.nv + p.bin.bal + p.bin.unbal + p.logNv
  beta.vec <- c(-beta, rep(beta, p))
  q <- p + 1
  error.fun <- function(e) return(list(coefficients = rep(NA, q)))
  ev <- tryCatch(eigen(crossprod(as.matrix(cbind(1, data[, -1]))), )$vectors[, q],
                 error = error.fun)
  angle <- as.numeric(acos(crossprod(ev, beta.vec) / 
                             (sqrt(crossprod(ev)) * sqrt(crossprod(beta.vec)))))
  epv <- sum(data$y) / p
  
  separation <- tryCatch(glm(y ~ ., family = "binomial", data = data,
                             method = "detect_separation"),
                         error = function(e) list(separation = NA))
  ml <- tryCatch(glm(y ~ ., family = "binomial", x = TRUE, data = data),
                 error = error.fun)
  cm <- tryCatch(glm(y ~ ., family = "binomial", data = data, method = "brglmFit",
                     type = "correction"),
                 error = error.fun)
  onest.jrle <- tryCatch(oneStepJRLE(y ~ ., data = data), error = error.fun)
  aurle <- tryCatch(oneStepJRLE(y ~ ., data = data, mle = ml), error = error.fun)
  firth <- tryCatch(glm(y ~ ., data = data, family = "binomial", method = "brglmFit",
                        type = "AS_mean"),
                    error = error.fun)
  kenne.pagui <- tryCatch(glm(y ~ ., data = data, family = "binomial",
                              method = "brglmFit",
                              type = "AS_median"), error = error.fun)
  ridge <- tryCatch(cv.glmnet(x = as.matrix(data[, colnames(data) != "y"]), y = data$y,
                              alpha = 0, family = "binomial"),
                    error = error.fun)
  if(p.nv == 0 && p.logNv == 0) {
    grouped.data <- groupData(data)
    grouped.data[, 1:2] <- grouped.data[, 1:2] + 0.5
    haldane <- tryCatch(glm(as.formula("cbind(Cases, Controls) ~ ."),
                            data = grouped.data, family = "binomial"), error = error.fun)
  } else {
    haldane <- list(coefficients = rep(NA, q), converged = NA)
  }
  estimates <- list(sep = separation$coefficients,
                    mle = ml$coefficients,
                    cm = cm$coefficients,
                    haldane = haldane$coefficients,
                    firth = firth$coefficients,
                    kennepagui = kenne.pagui$coefficients,
                    ridge = as.numeric(coefficients(ridge)),
                    oa = onest.jrle$coefficients,
                    wa = aurle$coefficients
                    )
  get.p.vals <- function(model) {
    tryCatch({
      p.vals <- summary(model)$coefficients[,"Pr(>|z|)"]
      # add NAs in case of non-convergence:
      aliased <- is.na(model$coefficients)
      if(any(aliased)) {
        tmp <- model$coefficients
        tmp[!aliased] <- p.vals
        tmp[aliased] <- NA
        p.vals <- tmp
      }
      return(p.vals)
    } , error = error.fun)
  }
  p.vals <- list(sep = rep(NA, q), 
                 mle = get.p.vals(ml),
                 cm = get.p.vals(cm),
                 haldane = get.p.vals(haldane),
                 firth = get.p.vals(firth),
                 kennepagui = get.p.vals(kenne.pagui),
                 ridge = rep(NA, q), 
                 oa = rep(NA, q),
                 wa = rep(NA, q))
  converged <- c(NA, ml$converged, cm$converged, haldane$converged, 
                 firth$converged, kenne.pagui$converged, NA, NA, NA)
  corr <- suppressWarnings(cor(cbind(1, data[, -1])))
  res <- data.frame(estimate = unlist(estimates), p.value = unlist(p.vals),
                    beta = beta.vec, 
                    epv = epv, angle = angle,
                    separation = separation$separation,
                    converged = rep(converged, each = q),
                    estimator = rep(c("Separation", "ML", "CM", "Haldane", "Firth",# "Exact", 
                                      "Kenne Pagui et al.", 
                                      "Ridge", "OA", "WA"), each = q),
                    trial = i, correlation  = corr)
  rownames(res) <- NULL
  return(res)
}