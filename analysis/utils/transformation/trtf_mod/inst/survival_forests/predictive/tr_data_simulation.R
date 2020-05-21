library("survival")
nsim <- 100
ntrain <- 250
ntest <- 500
set.seed(33)

### scale and shape parameters of Weibull distribution
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

### data-dependent shift parameter of the model
shift_f <- function(d, prod_shift = c(1, 1), ...){
  with(d, f1(X1, X2, X3, X4, X5, scale = TRUE)) * prod_shift[[1]] + 
          with(d, f1(X6, X7, X8, X9, X10, scale = TRUE)) * prod_shift[[2]] * d$trt}

### data-dependent scale parameter of the model
scale_f <- function(d, prod_scale = c(1, 1), ...){
  exp(with(d, f1(X11, X12, X13, X14, X15, scale = TRUE)) * prod_scale[[1]] + 
      with(d, f1(X16, X17, X18, X19, X20, scale = TRUE)) * prod_scale[[2]] * d$trt)}

### Model of the outcome
rmodel <- function(n = 100, shift = 0, scale = 1) {
  ### Gompertz distribution
  FZ <- mlt:::.MinExtrVal()   
  q0 <- function(u) qweibull(u, scale = wscale, shape = wshape)
  q1 <- function(u, shift, scale) q0(FZ$p((FZ$q(u) - shift) / scale))
  y <- q1(runif(n), shift = shift, scale = scale)

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

ry <- function(d, prod_shift, prod_scale,...)
  rmodel(nrow(d),
         shift = shift_f(d, prod_shift = prod_shift), 
         scale = scale_f(d, prod_scale = prod_scale))


### data generation
dgp <- function(ntrain, ntest, p = 5, prod_shift = c(0, 0), prod_scale = c(0, 0), ...) {
  n <- ntrain + ntest
  d <- data.frame(matrix(runif(20 * n), ncol = 20),
                  matrix(runif(n * p, min = 0, max = 1), nrow = n))
  d$trt <- sample(rep(0:1, length = n))
  d$y <- ry(d, prod_shift=prod_shift, prod_scale=prod_scale)
  d$y <- with(d, Surv(y, rep(1, nrow(d))))
  list('train' = d[1:ntrain,], 'test' = d[(ntrain+1):n,],
       prod_shift = prod_shift, prod_scale = prod_scale)
}

### True model log-likelihood
mylogLik <- function(d){
  shift <- shift_f(d$test, prod_shift = d$prod_shift)
  scale <- scale_f(d$test, prod_scale = d$prod_scale)
  pw <- function(u) pweibull(u, scale=wscale, shape=wshape)
  FZ <- mlt:::.MinExtrVal()
  h <- function(u) FZ$q(pw(u))
  dh <- function(u) - dweibull(u, scale=wscale, shape=wshape) / ((1 - pw(u)) * log(1 - pw(u)))
  ty <- h(d$test$y[,1]) * scale + shift
  dhy <- dh(d$test$y[,1])
  dh_idx <- !is.na(dhy) 
  sum(ty[dh_idx] - exp(ty[dh_idx]) + log(scale[dh_idx]) + log(dhy[dh_idx]))
}

if (FALSE) {
### some basic checks
library("tram")
N <- 10000

d <- dgp(N, N, prod_shift = c(0, 0), prod_scale = c(0, 0))
m0 <- as.mlt(Survreg(y ~ 1, data = d$train))
logLik(m0, newdata = d$test)
mylogLik(d)
coef(m0)

d <- dgp(N, N, prod_shift = c(1, 0), prod_scale = c(0, 0))
d$train$sh <- shift_f(d$train, prod_shift = d$prod_shift)
d$test$sh <- shift_f(d$test, prod_shift = d$prod_shift)
m1 <- as.mlt(Survreg(y ~ sh, data = d$train, fixed = c("sh" = -1)))
logLik(m1, newdata = d$test, parm = coef(m1))
mylogLik(d)
coef(m1)

d <- dgp(N, N, prod_shift = c(0, 1), prod_scale = c(0, 0))
d$train$sh <- shift_f(d$train, prod_shift = d$prod_shift)
d$test$sh <- shift_f(d$test, prod_shift = d$prod_shift)
m2 <- as.mlt(Survreg(y ~ sh, data = d$train, fixed = c("sh" = -1)))
logLik(m2, newdata = d$test, parm = coef(m2))
mylogLik(d)
coef(m2)

d <- dgp(N, N, prod_shift = c(1, 1), prod_scale = c(0, 0))                                                                                                                                                  
d$train$sh <- shift_f(d$train, prod_shift = d$prod_shift)
d$test$sh <- shift_f(d$test, prod_shift = d$prod_shift)
m3 <- as.mlt(Survreg(y ~ sh, data = d$train, fixed = c("sh" = -1)))
logLik(m3, newdata = d$test, parm = coef(m3))
mylogLik(d)
coef(m3)

d <- dgp(N, N, prod_shift = c(0, 0), prod_scale = c(1, 0))                                                                                                                                                 
d$train$sh <- shift_f(d$train, prod_shift = d$prod_shift)
d$train$sc <- scale_f(d$train, prod_scale = d$prod_scale)
d$test$sh <- shift_f(d$test, prod_shift = d$prod_shift)
d$test$sc <- scale_f(d$test, prod_scale = d$prod_scale)
m4 <- as.mlt(Survreg(y | sc + 0 ~ sh, data = d$train, fixed = c("sh" = 0)))                                                                                                                                 
logLik(m4, newdata = d$test, parm = coef(m4))
mylogLik(d)
coef(m4)

d <- dgp(N, N, prod_shift = c(0, 0), prod_scale = c(0, 1))                                                                                                                                                  
d$train$sh <- shift_f(d$train, prod_shift = d$prod_shift)
d$train$sc <- scale_f(d$train, prod_scale = d$prod_scale)
d$test$sh <- shift_f(d$test, prod_shift = d$prod_shift)
d$test$sc <- scale_f(d$test, prod_scale = d$prod_scale)
m5 <- as.mlt(Survreg(y | sc + 0 ~ sh, data = d$train, fixed = c("sh" = 0)))                                                                                                                                 
logLik(m5, newdata = d$test, parm = coef(m5))                                                                                                                                                               
mylogLik(d)
coef(m5)

d <- dgp(N, N, prod_shift = c(0, 0), prod_scale = c(1, 1))
d$train$sh <- shift_f(d$train, prod_shift = d$prod_shift)
d$train$sc <- scale_f(d$train, prod_scale = d$prod_scale)
d$test$sh <- shift_f(d$test, prod_shift = d$prod_shift)
d$test$sc <- scale_f(d$test, prod_scale = d$prod_scale)
m6 <- as.mlt(Survreg(y | sc + 0 ~ sh, data = d$train, fixed = c("sh" = 0)))
logLik(m6, newdata = d$test, parm = coef(m6))
mylogLik(d)
coef(m6)

d <- dgp(N, N, prod_shift = c(1, 1), prod_scale = c(1, 1))
d$train$sh <- shift_f(d$train, prod_shift = d$prod_shift)
d$train$sc <- scale_f(d$train, prod_scale = d$prod_scale)
d$test$sh <- shift_f(d$test, prod_shift = d$prod_shift)
d$test$sc <- scale_f(d$test, prod_scale = d$prod_scale)
m7 <- as.mlt(Survreg(y | sc + 0 ~ sh, data = d$train, fixed = c("sh" = -1)))
logLik(m7, newdata = d$test, parm = coef(m7))
mylogLik(d)
coef(m7)

}

args <- data.frame(p = rep(c(5, 50), times = 3),
                   prod_shift_int = rep(c(1, 0, 1), each = 2),
                   prod_shift_trt = rep(c(1, 0, 1), each = 2),
                   prod_scale_int = rep(c(0, 1, 1), each = 2),
                   prod_scale_trt = rep(c(0, 1, 1), each = 2))

simdat <- vector(mode = "list", length = nrow(args))
loglik <- vector(mode = "list", length = nrow(args))

for (i in 1:nrow(args)) {

    simdat[[i]] <- replicate(nsim, dgp(ntrain = ntrain,
                                       ntest = ntest,
                                       p = args[i, "p"],
                                       prod_shift = c(args[i, "prod_shift_int"], args[i, "prod_shift_trt"]),
                                       prod_scale = c(args[i, "prod_scale_int"], args[i, "prod_scale_trt"])),
                             simplify = FALSE)
    loglik[[i]] <- lapply(simdat[[i]], mylogLik)
}

save(args, simdat, loglik, file='tr_simulated_data.rda')
