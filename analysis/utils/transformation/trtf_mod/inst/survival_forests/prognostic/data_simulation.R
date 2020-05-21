
library("survival")
## number of repetitions
nsim <- 100
## size of learning sample
ntrain <- 250
## size of validation sample
ntest <- 500
set.seed(33)

## scale and shape parameters of Weibull distribution
wscale <- 1
wshape <- 1

### Friedman function
f1 <- function(x1, x2, x3, x4, x5, scale = TRUE) {
  ret <- 10 * sin(pi * x1 * x2) + 20 * (x3 - 0.5)^2 + 10 * x4 + 5 * x5
  if (scale) {
    ret <- ret - min(ret) 
    ret <- ret / max(ret)
    ret <- 3 * ret - 1.5
  }
  ret
}

### model scale function
## if there is no scale in the model, prod_scale = 0
scale_f <- function(d, prod_scale = 1, ...){
  exp(with(d, f1(X1, X2, X3, X4, X5, scale = TRUE)) * prod_scale)}

### model shift function
## if there is no shift in the model, prod_shift = 0
shift_f <- function(d, prod_shift = 1, ...){
  with(d, f1(X6, X7, X8, X9, X10, scale = TRUE)) * prod_shift}

rmodel <- function(n = 100, shift = 0, scale = 1) {
  ### Gompertz distribution
  FZ <- mlt:::.MinExtrVal()   
  q0 <- function(u) qweibull(u, scale = wscale, shape = wshape)
  q1 <- function(u, shift, scale) q0(FZ$p((FZ$q(u) - shift) / scale))
  y <- q1(runif(n), shift = shift, scale = scale)
  
  ### to fix infinte generated values
  inf_idx <- is.infinite(y)
  ninf <- sum(inf_idx)
  while (ninf > 0){
    y[inf_idx] <- q1(runif(ninf), shift = shift[inf_idx], 
                     scale = scale[inf_idx]) 
    inf_idx <- is.infinite(y) 
    ninf <- sum(inf_idx)
  }
  return(pmax(.Machine$double.eps, y))
}

ry <- function(d, prod_shift, prod_scale, ...)
  rmodel(nrow(d), shift_f(d, prod_shift=prod_shift), 
         scale_f(d, prod_scale=prod_scale))


### data generation
dgp <- function(ntrain, ntest, p = 5, prod_shift = 0, prod_scale = 0, ...) {
  n <- ntrain + ntest
  d <- data.frame(matrix(runif(10 * n), ncol = 10),
                  matrix(runif(n * p, min = 0, max = 1), nrow = n))
  d$y <- ry(d, prod_shift=prod_shift, prod_scale=prod_scale)
  d$y <- with(d, Surv(y, rep(1, nrow(d))))
  list('train' = d[1:ntrain,], 'test' = d[(ntrain+1):n,],
       prod_shift = prod_shift, prod_scale = prod_scale)
}

### log-likelihood of the model
mylogLik <- function(d){
  shift <- shift_f(d$test, prod_shift = d$prod_shift)
  scale <- scale_f(d$test, prod_scale = d$prod_scale)
  pw <- function(u) pweibull(u, scale = wscale, shape = wshape)
  FZ <- mlt:::.MinExtrVal()
  h <- function(u) FZ$q(pw(u))
  dh <- function(u) - dweibull(u, scale = wscale, shape = wshape) / ((1 - pw(u)) * log(1 - pw(u)))
  ty <- scale * h(d$test$y[,1]) + shift 
  sum(ty - exp(ty) + log(scale) + log(dh(d$test$y[,1])))
}

if (FALSE) {
### some basic checks
library("tram")
N <- 50000

d <- dgp(N, N, prod_shift = 0, prod_scale = 0)
m0 <- as.mlt(Survreg(y ~ 1, data = d$train))
logLik(m0, newdata = d$test, parm = coef(m0))
logLik(m0, newdata = d$test, parm = c(0, 1))
mylogLik(d)
coef(m0)

d <- dgp(N, N, prod_shift = 1, prod_scale = 0)
d$train$sh <- shift_f(d$train)
d$test$sh <- shift_f(d$test)
m1 <- as.mlt(Survreg(y ~ sh, data = d$train, fixed = c("sh" = -1)))
logLik(m1, newdata = d$test, parm = coef(m1))
logLik(m1, newdata = d$test, parm = c(0, 1, -1))
mylogLik(d)
coef(m1)

d <- dgp(N, N, prod_shift = 1, prod_scale = 1)
d$train$sh <- shift_f(d$train)
d$train$sc <- scale_f(d$train)
d$test$sh <- shift_f(d$test)
d$test$sc <- scale_f(d$test)
d$train$int <- 1
d$test$int <- 1
m2 <- as.mlt(Survreg(y | sc + 0 ~ int + sh, data = d$train, 
                     fixed = c("(Intercept):sc" = 0, "sh" = -1)))
logLik(m2, newdata = d$test, parm = coef(m2))
logLik(m2, newdata = d$test, parm = c(0, 1, 0, -1))
mylogLik(d)
coef(m2)

d <- dgp(N, N, prod_shift = 0, prod_scale = 1)
d$train$sh <- shift_f(d$train)
d$train$sc <- scale_f(d$train)
d$test$sh <- shift_f(d$test)
d$test$sc <- scale_f(d$test)
d$train$int <- 1
d$test$int <- 1
m3 <- as.mlt(Survreg(y | sc + 0 ~ int + sh, data = d$train, 
                     fixed = c("(Intercept):sc" = 0, "sh" = 0)))
logLik(m3, newdata = d$test, parm = coef(m3))
logLik(m3, newdata = d$test, parm = c(0, 1, 0, 0))
mylogLik(d)
coef(m3)

}

args <- expand.grid(p = c(5, 50), 
                    prod_shift = c(0, 1), 
                    prod_scale = c(0, 1))

simdat <- vector(mode = "list", length = nrow(args))
loglik <- vector(mode = "list", length = nrow(args))

for (i in 1:nrow(args)) {

    simdat[[i]] <- replicate(nsim, dgp(ntrain = ntrain,
                                       ntest = ntest,
                                       p = args[i, "p"],
                                       prod_shift = args[i, "prod_shift"],
                                       prod_scale = args[i, "prod_scale"]),
                             simplify = FALSE)
    loglik[[i]] <- lapply(simdat[[i]], mylogLik)
}

save(args, simdat, loglik, file='simulated_data.rda')
