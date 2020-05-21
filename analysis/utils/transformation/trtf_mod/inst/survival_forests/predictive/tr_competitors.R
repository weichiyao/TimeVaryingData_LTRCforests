library("trtf")
library("tram")
library("party")
library("parallel")
library("survival")
set.seed(7782)

### method for for trtf::traforest
predict.mytrtf <- function(object, newdata = object$data, OOB = FALSE, type = "weights", ...)
{
    stopifnot(!OOB)
    predict(object$trtf, newdata = newdata, OOB = FALSE, type = "weights")
}

### joint parameters for all forest procecdures
## number of trees
ntree <- 250

## maximal depth of trees
nodedepth <- 10
## number of obs in terminal node
minbucket <- 20
## number of variables to choose from in a node

### parameters of Bernstein approximation
## polynom order
ord <- 5
## whether approximation applied to the log-transformed time or not
logf <- TRUE

### parameters for party/partykit
ctrl_partykit <- partykit:::ctree_control(teststat = "quad",
                                          testtype = "Univ",
                                          mincriterion = 0,
                                          minbucket = minbucket,
                                          maxdepth = nodedepth)


### Weibull distributional survival forest 
### with log-rank prognosic and predictive splitting
### W(alpha, beta)
w_alpha_beta <- function(trtf, learn, test, strata){
    mtry <- trtf$mtry
    trtf <- trtf$object

    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"
    
    m0 <- Survreg(y ~ trt, data = learn)
    rf <- traforest(as.mlt(m0),
                    formula = y | trt ~ .,
                    data = learn,
                    parm = c("(Intercept)", "trt"),
                    mtry = mtry,
                    ntree = ntree,
                    strata = strata,
                    control = ctrl_partykit,
                    weights = initw)

    ### use hack to wire rf weights into predict.traforest
    trtf$trtf <- rf
    class(trtf) <- c(class(trtf)[1], "mytrtf", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(loglik = llNN, object = rf)
}

### Weibull distributional survival forest
### with general prognostic splitting and log-rank predictive splitting
### W(theta, beta)
w_theta_beta <- function(trtf, learn, test, strata){
    mtry <- trtf$mtry
    trtf <- trtf$object
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"

    m0 <- Survreg(y ~ trt, data = learn)
    rf <- traforest(as.mlt(m0),
                    formula = y | trt ~ .,
                    data = learn,
                    mtry = mtry,
                    ntree = ntree,
                    strata = strata,
                    control = ctrl_partykit,
                    weights = initw)
    ### use hack to wire rf weights into predict.traforest
    trtf$trtf <- rf
    class(trtf) <- c(class(trtf)[1], "mytrtf", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(loglik = llNN, object = rf)
}

### Weibull distributional survival forest
### with general prognostic and predictive splitting
### W(theta, theta)
w_theta_theta <- function(trtf, learn, test, strata){
    mtry <- trtf$mtry
    trtf <- trtf$object
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"

    m0 <- Survreg(y | trt ~ 1, data = learn)
    rf <- traforest(as.mlt(m0),
                    formula = y | trt ~ .,
                    data = learn,
                    mtry = mtry,
                    ntree = ntree,
                    strata = strata,
                    control = ctrl_partykit,
                    weights = initw)
    ### use hack to wire rf weights into predict.traforest
    trtf$trtf <- rf
    class(trtf) <- c(class(trtf)[1], "mytrtf", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(loglik = llNN, object = rf)
}

### Transformation survival forest with Bernstein approximation and
### log-rank prognostic and predictive splitting
### Bs(alpha, beta)
bs_alpha_beta <- function(trtf, learn, test, strata, order = ord, log_first = logf){

    mtry <- trtf$mtry
    trtf <- trtf$object
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"

    ## since bernstein basis does not include an intercept
    ## we add it to the data implicitly
    learn$int <- 1
    test$int <- 1

    ## the largest time values are not used in approximation
    sup <- c(sqrt(.Machine$double.eps), quantile(learn$y[,1], .95))
    m0 <- Coxph(y ~ int + trt, 
                data = learn, 
                fixed = c("int"), 
                order = order, 
                log_first = log_first, 
                support = sup)
    rf <- traforest(as.mlt(m0),
                    formula = y | int + trt ~ .,
                    data = learn,
                    parm = c("int", "trt"),
                    mltargs = list(fixed = c("int" = 0)),
                    mtry = mtry,
                    ntree = ntree,
                    strata = strata,
                    control = ctrl_partykit,
                    weigths = initw)
    ### use hack to wire rf weights into predict.traforest
    trtf$trtf <- rf
    class(trtf) <- c(class(trtf)[1], "mytrtf", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(loglik = llNN, object = rf)
}

### Transformation survival forest with Bernstein approximation
### wtih general prognostic and log-rank predictive splitting
### BS(theta, beta)
bs_theta_beta <- function(trtf, learn, test, strata, order = ord, log_first = logf){
    mtry <- trtf$mtry
    trtf <- trtf$object
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"

    ## the largest time values are not used in approximation                  
    sup <- c(sqrt(.Machine$double.eps), quantile(learn$y[,1], .95))
    m0 <- Coxph(y ~ trt, 
                data = learn, 
                order = order, 
                log_first = log_first, 
                support = sup)
    rf <- traforest(as.mlt(m0),
                    formula = y | trt ~ .,
                    data = learn,
                    mtry = mtry,
                    ntree = ntree,
                    strata = strata,
                    control = ctrl_partykit,
                    weights = initw)
    ### use hack to wire rf weights into predict.traforest
    trtf$trtf <- rf
    class(trtf) <- c(class(trtf)[1], "mytrtf", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(loglik = llNN, object = rf)
}

### Transformation survival forest with Bernstein approximation
### with general prognostic and predictive splitting
### BS(theta, theta)
bs_theta_theta <- function(learn, test, strata, mtry, order = ord, log_first = logf){
    ## we do not pass any learned forest to this function since we use it
    ## for learning the first forest which is used for learning other forests
    ## the largest time values are not used in approximation                  
    sup <- c(sqrt(.Machine$double.eps), quantile(learn$y[,1], .95))
    m0 <- Coxph(y | trt ~ 1, 
                data = learn, 
                order = order, 
                log_first = log_first, 
                support = sup)
    rf <- traforest(as.mlt(m0),
                    formula = y | trt ~ .,
                    data = learn,
                    mtry = mtry,
                    ntree = ntree,
                    strata = strata,
                    control = ctrl_partykit)
    ### nearest neighbor predictions
    cf <- predict(rf, newdata = test, type = "coef")
    llNN <- logLik(rf, newdata = test, coef = cf)

    list(loglik = llNN, object = rf, mtry = mtry)
}

### Cox model with treatment and treatment-variable interactions
tcoxph <- function(learn, test, order = ord, log_first = logf){
    ## make model with treatment and treatment-variable interactions
    fm <- names(learn)
    fm <- paste(fm[!(fm %in% c("y", "trt"))], collapse = "+")
    fmy <- as.formula(paste("y ~ trt +", fm, " + trt * (", fm, ")"))

    ## the largest time values are not used in approximation
    sup <- c(sqrt(.Machine$double.eps), quantile(learn$y[,1], .95))
    cph <- Coxph(fmy , learn, order = order, log_first = log_first, support = sup)
    list(loglik = logLik(cph, newdata = test), object = cph)
}

### Weibull model with treatment and treatment-varibale interactions
tweib <- function(learn, test){
    ## make model with treatment and treatment-variable interactions
    fm <- names(learn)
    fm <- paste(fm[!(fm %in% c("y", "trt"))], collapse = "+")
    fmy <- as.formula(paste("y ~ trt +", fm, " + trt * (", fm, ")"))

    tw <- Survreg(fmy, learn)
    list(loglik = logLik(tw, newdata = test), object = tw)
}
