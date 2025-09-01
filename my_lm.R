# a simple helper function
matrix_trace <- function(mtx) {
  return(sum(diag(mtx)))
}

my_lm <- function(y, X, alpha=0.05) {
  
  if (is.vector(y)) {
    y <- matrix(y, ncol=1)
  }
  
  # set up the atomic pieces that we're going to use multiple times 
  y_t <- t(y)
  X_t <- t(X)
  scaling <- solve(X_t%*%X)
  n <- nrow(y)
  
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
  SSE <- y_t%*%y - y_t%*%H_f%*%y
  df_e <- n - matrix_trace(H_f)
  MSE <- (SSE/df_e) |> as.numeric()
  
  # estimate the covariance matrix of the parameter estimates
  covB <- MSE*scaling
  # standard errors of the parameter estimates
  seB <- sqrt(diag(covB))
  
  quad_residuals <- (H_f - H_r)
  
  SSR <- y_t%*%quad_residuals%*%y
  MSR <- SSR/matrix_trace(quad_residuals)
  
  tests <- test.o(B, seB, alpha, df_e)
  
  df <- c(
    n - 1,
    df_e,
    n - 1 - df_e
  )
  
  mtx_estimates <- cbind(
    "estimates" = B,
    "stderr" = seB,
    "t" = tests$t,
    "p" = tests$pvals,
    "CI_lwr" = tests$lwr,
    "CI_upr" = tests$upr
  )
  
  colnames(mtx_estimates) <- c("Estimate", "Std err", "t", "Pr(>|t|)", "2.5%", "97.5%")
  
  return(
    list(
      "coefficients" = B,
      "summary_table" = mtx_estimates,
      "residuals" = resids,
      "cov" = covB,
      "MSE" = MSE,
      "MSR" = MSR,
      "SSTO" = SSR + SSE,
      "df" = df
    )
  )
}

my_anova <- function(y, X_f, X_r=NULL) {
  if(is.vector(y)) {
    y <- matrix(y, ncol=1)
  }
  n <- nrow(y)
  
  # create the transposes that we'll need
  y_t <- t(y)
  X_f_t <- t(X_f)
  X_r_t <- t(X_r)
  
  # compute hat matrices for the full and reduced models
  H_f <- X_f%*%solve(X_f_t%*%X_f)%*%X_f_t
  H_r <- X_r%*%solve(X_r_t%*%X_r)%*%X_r_t
  
  # now compute the Mean Deviation of the two models and the MSE of the full
  # model
  df <- c(matrix_trace(H_f - H_r), n - matrix_trace(H_f))
  SSD <- (y_t%*%(H_f - H_r)%*%y)
  MSD <- SSD/df[1]
  SSE_f <- y_t%*%(y - H_f%*%y)
  MSE_f <- SSE_f/df[2]
  SSTO <- sum((y-mean(y))^2)
  
  return(
    list(
      "F" = MSD/MSE_f,
      "df" = df,
      "MSD" = MSD,
      "SSE_F" = SSE_f,
      "SSD" = SSD,
      "SSTO" = SSTO
    )
  )
}

# perform hypothesis testing on the estimates
test.o <- function(B, seB, alpha, df_e) {
  t <- B/seB
  pvals <- 2*pt(abs(t), df=df_e, lower.tail=FALSE)
  
  t_crit <- qt(1-alpha/2, df=df_e)
  
  CI_lwr <- B - t_crit*seB
  CI_upr <- B + t_crit*seB
  
  return(
    list(
      "t" = t,
      "lwr" = CI_lwr,
      "upr" = CI_upr,
      "pvals" = pvals
    )
  )
}
