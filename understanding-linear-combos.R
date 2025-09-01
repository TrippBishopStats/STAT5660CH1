rm(list=ls())
library(matlib)
source("my_lm.R")

df_cereal <- read_csv("KentonFoodData.csv") |> 
  mutate(
    store = as_factor(store),
    design = as_factor(design)
  )

df_cereal <- df_cereal |> na.omit()

y <- df_cereal$sales
X_f <- model.matrix(~ design, data=df_cereal)

my.fit <- my_lm(y, X_f)

l <- matrix(c(0,-1,1,0), ncol=1)
cov_B <- my.fit$cov
B <- my.fit$coefficients

# this is the same thing y3_bar - y2_bar. These are (y3_bar - y1_bar) - (y2_bar
# - y1_bar) specifically.
t(l)%*%B

# let's look at the std errs
t(l)%*%cov_B%*%l |> sqrt()

t(l)%*%B/(t(l)%*%cov_B%*%l |> sqrt())


t(l)%*%cov_B[,1, drop=FALSE]
t(l)%*%cov_B[,2, drop=FALSE]
t(l)%*%cov_B[,3, drop=FALSE]
t(l)%*%cov_B[,4, drop=FALSE]
