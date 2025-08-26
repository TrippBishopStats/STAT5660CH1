rm(list=ls())
y <- c(85,95,69,58,41,74,71,52,41,34,50,40)
xu <- c(60,72,61,50,54,68,66,59,56,56,55,51)
# centered data
x <- xu-mean(xu) # use x! #
t <- as.factor((1:3)%x%rep(1,4)) # rep(1:3, each=4) does the same thing
n <- length(y)

lm(y ~ t + x)

# t1 is the referent group
X <- model.matrix(~ t + x)