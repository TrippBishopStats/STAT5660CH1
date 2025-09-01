rm(list=ls())
library(tidyverse)
source("my_lm.R")

df_cereal <- read_csv("KentonFoodData.csv") |> 
  mutate(
    store = as_factor(store),
    design = as_factor(design)
  )

df_cereal <- df_cereal |> na.omit()

lm.fit <- lm(sales~design, data=df_cereal)

y <- df_cereal$sales
X_f <- model.matrix(~ design, data=df_cereal)
X_r <- model.matrix(~ 1, data=df_cereal)

my.fit <- my_lm(y, X_f)

### ANOVA output
my_anova(y, X_f, X_r)
anova(lm.fit)

1+1










l <- matrix(c(0,-1,1,0), ncol=1)
cov_B <- my.fit$cov

est <- my.fit$estimates[,1, drop=FALSE]

(t(l)%*%est)/(t(l)%*%cov_B%*%l |> sqrt())

l <- matrix(c(1,-1,0,0,
              1,0,-1,0,
              1,0,0,-1), nrow=3, ncol=4, byrow=TRUE)

l <- matrix(c(0,-1,1,0), nrow=1, ncol=4, byrow=TRUE)

library(multcomp)
c.fit <- glht(lm.fit, linfct=l)
c.sum <- summary(c.fit, test=univariate())$test
c.t <- cbind(c.sum$coefficients, c.sum$sigma, c.sum$tstat, c.sum$pvalues)
c.ci <- confint(c.fit, calpha=univariate_calpha())$confint[,-1]

# compute the estimates of the contrasts
l%*%est

# compute the standard errors
SE <- diag(l%*%cov_B%*%t(l))|> sqrt()

# these two lines will compute the t statistics

D <- SE |> diag(ncol=3,nrow=3) |> solve()
t(l%*%est)%*%D
              