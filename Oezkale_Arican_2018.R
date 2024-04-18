oneStepJRLE <- function(formula = NULL, data = NULL, start = NULL, 
                        control = list(maxit = 1),
                        mle = glm(formula, data = data, start = start, 
                                  control = control, family = binomial(),
                                  x = TRUE)){
  X <- mle$x
  beta <- mle$coefficients
  p <- length(beta)
  W <- diag(mle$fitted.values * (1 - mle$fitted.values))
  mat1 <- t(X) %*% W %*% X
  eig <- eigen(mat1)
  alphas <- (t(eig$vectors) %*% beta)
  kis <- (1 + sqrt(1 + eig$values * alphas^2)) / alphas^2 
  k <- p / sum(1/kis)
  mat <- solve(mat1 + k * diag(p))
  beta.hat <- drop((diag(p) - k^2 * mat %*% mat) %*% beta)
  return(list(coefficients = beta.hat))
}
