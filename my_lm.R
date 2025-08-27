my_lm <- function(y, X, alpha=0.05) {
  if (!is.vector(y)) {
    stop("y must be a vector")
  } else if (!is.matrix(X)) {
    stop("X must be a matrix")
  } else if(alpha <= 0 | alpha >= 1) {
    stop("confidence level must be between 0 and 1")
  }
  
  # set up the atomic pieces that we're going to use multiple times 
  y <- matrix(y, ncol=1) # make the response a proper column vector
  y_t <- t(y)
  X_t <- t(X)
  scaling <- solve(X_t%*%X)
  n <- nrow(y)
  I <- diag(1, ncol=n, nrow=n)
  
  # the reduced model is just a mean (null) model so the hat matrix is
  # equivalent to 1/n*J
  H_r <- 1/n*matrix(rep(1, times=n*n), ncol=n, nrow=n)
  
  # estimate the population parameters
  B <- scaling%*%X_t%*%y
  H_f <- X%*%scaling%*%X_t
  
  # compute the residuals
  resids <- y - H_f%*%y
  
  # we need to estimate the mean square error first, so that we can use the
  # result in other estimates.
  quad_err <- (I - H_f)
  
  SSE <- y_t%*%quad_err%*%y
  MSE <- (SSE/sum(diag(quad_err))) |> as.numeric()
  
  # estimate the covariance matrix of the parameter estimates
  covB <- MSE*scaling
  # standard errors of the parameter estimates
  seB <- sqrt(diag(covB))
  
  quad_residuals <- (H_f - H_r)
  
  SSR <- y_t%*%quad_residuals%*%y
  MSR <- SSR/sum(diag(quad_residuals))
  
  t <- B/seB
  pvals <- 2*pt(abs(t), df=(n-1), lower.tail=FALSE)
  
  t_crit <- qt(1-alpha/2, df=(n-1))
  
  CI_lwr <- B - t_crit*seB
  CI_upr <- B + t_crit*seB
  
  mtx_estimates <- cbind(
    "estimates" = B,
    "stderr" = seB,
    "t" = t,
    "p" = pvals,
    "CI_lwr" = CI_lwr,
    "CI_upr" = CI_upr 
  )
  
  colnames(mtx_estimates) <- c("Estimate", "Std err", "t", "Pr(>|t|)", "2.5%", "97.5%")
  
  return(
    list(
      "estimates" = mtx_estimates,
      "residuals" = resids,
      "MSE" = MSE,
      "MSR" = MSR
    )
  )
}

# test harness
y <- df_cereal$sales

my_lm(y, X)

# lm(sales~design, data=df_cereal) |> summary()
