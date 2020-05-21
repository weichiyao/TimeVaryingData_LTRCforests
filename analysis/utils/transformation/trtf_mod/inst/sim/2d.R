
library("trtf")

nsim <- 100
ntree <- 100

sigma <- function(d, tau = .01, prod_sigma = 1, ...)
    1 + pnorm(d$x1, sd = tau, mean = .5) * prod_sigma

mu <- function(d, tau = .01, prod_mu = 1, ...)
    pnorm(d$x2, sd = tau, mean = .5) * prod_mu

med <- function(d, ...)
    mu(d, ...)

qu <- function(d, prob, prod_mu, prod_sigma, tau, ...)
    do.call("cbind", lapply(prob, function(p) qnorm(p, mean = mu(d, tau = tau, prod_mu = prod_mu), 
                                                    sd = sigma(d, tau = tau, prod_sigma = prod_sigma))))

ry <- function(d, tau, prod_mu, prod_sigma, ...)
    rnorm(nrow(d), mean = mu(d, tau = tau, prod_mu = prod_mu), sd = sigma(d, tau = tau, prod_sigma = prod_sigma))

dgp <- function(n = 250, p = 5, ...) {
    d <- data.frame(x1 = runif(n), 
                    x2 = runif(n),
                    matrix(runif(n * p), nrow = n))
    d$y <- ry(d, ...)
    d	
}

loglik <- function(d, p, tau, prod_mu, prod_sigma, ...)
    sum(dnorm(d$y, mean = mu(d, tau = tau, prod_mu = prod_mu), sd = sigma(d, tau = tau, prod_sigma = prod_sigma), log = TRUE))

parm <- expand.grid(nsim = nsim, n = 250, ntest = 250, p = c(5, 500), tau = 0.01, 
                    prod_sigma = c(0, 1), prod_mu = c(0, 1))
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
#                                            "tfbagg", 
                                            "tfrf", 
#                                            "tfrfExSplit" 
#                                            "tfbaggBern", 
                                            "tfrfBern", 
#                                            "cfbagg", 
                                            "cfrf", 
#                                            "bagg", 
                                            "rf",
                                            "rfBern"
))
out <- mclapply(parm, function(x) do.call("sf", x), mc.cores = length(parm))

save(out, file = "2d.rda")
