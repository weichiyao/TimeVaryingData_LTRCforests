
### same as friedman.R

library("trtf")

nsim <- 10
ntree <- 100

f1 <- function(x1, x2, x3, x4, x5, scale = TRUE) {
    ret <- 10 * sin(pi * x1 * x2) + 20 * (x3 - 0.5)^2 + 10 * x4 + 5 * x5
    if (scale) {
        ret <- ret - min(ret) 
        ret <- ret / max(ret)
        ret <- 3 * ret - 1.5
    }
    ret
}

sigma <- function(d, prod_sigma = 1, ...)
    exp(with(d, f1(X1, X2, X3, X4, X5, scale = TRUE)) * prod_sigma)

mu <- function(d, prod_mu = 1, ...)
    with(d, f1(X6, X7, X8, X9, X10, scale = TRUE)) * prod_mu

med <- function(d, ...)
    mu(d, ...)

qu <- function(d, prob, prod_mu, prod_sigma, ...)
    do.call("cbind", lapply(prob, function(p) qnorm(p, mean = mu(d, prod_mu = prod_mu), 
                                                    sd = sigma(d, prod_sigma = prod_sigma))))

ry <- function(d, prod_mu, prod_sigma, ...)
    rnorm(nrow(d), mean = mu(d, prod_mu = prod_mu), sd = sigma(d, prod_sigma = prod_sigma))

dgp <- function(n = 250, p = 5, ...) {
    d <- data.frame(matrix(runif(10 * n), ncol = 10),
                    matrix(runif(n * p, min = 0, max = 1), nrow = n))
    d$y <- ry(d, ...)
    d	
}


loglik <- function(d, p, prod_mu, prod_sigma, ...)
    sum(dnorm(d$y, mean = mu(d, prod_mu = prod_mu), sd = sigma(d, prod_sigma = prod_sigma), log = TRUE))

parm <- expand.grid(nsim = nsim, n = c(250, 500, 1000, 2000, 4000), 
                    ntest = 10, p = 5, tau = 0, prod_sigma = 1, prod_mu = 1)
parm <- lapply(1:nrow(parm), function(i) parm[i,])

source("artificial_simfun.R")

ctrl <- ctree_control(teststat = "quad", testtype = "Univ", 
                      mincriterion = 0, saveinfo = FALSE, 
                      minsplit = 25, minbucket = 10)

mc.reset.stream()
sf <- function(...) simfun(..., methods = c("ct", 
                                            "tt", 
                                            "ttExSplit",
                                            "ttM", 
                                            "ttBern", 
                                            "ttBernExSplit",
                                            "tfbagg", 
                                            "tfrf", 
#                                            "tfrfExSplit" 
                                            "tfbaggBern", 
                                            "tfrfBern", 
                                            "cfbagg", 
                                            "cfrf", 
                                            "bagg", 
                                            "rf"
))

parm

out <- mclapply(parm, function(x) do.call("sf", x), mc.cores = length(parm))

save(out, file = "timings.rda")
