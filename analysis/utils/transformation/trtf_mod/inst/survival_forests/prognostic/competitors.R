
library("trtf")
library("party")
library("tram")
library("survival")
set.seed(29)

### joint parameters for all forest procecdures
## number of trees
ntree <- 250
## maximal depth of trees
nodedepth <- 10
## number of obs in terminal node
minbucket <- 20

### parameters for Bernstein-based forest procedures
## order of approximation
ord <- 5
## whather to make log-transformation before the approximation
logf <- TRUE

### parameters for party/partykit
ctrl_partykit <- partykit:::ctree_control(teststat = "quad", 
                                          testtype = "Univ", 
                                          mincriterion = 0,
                                          minbucket = minbucket,
                                          maxdepth = nodedepth)

### these are _HACKS_ to extract NN weights from various forests
### trtf::predict.traforest calls these methods to obtain the nearest
### neighbor weights and then proceeds to estimate the coefficients
## method for party::cforest
predict.mycf <- function(object, newdata = object$data, OOB = FALSE, type = "weights", ...)
{
    stopifnot(!OOB)
    do.call("cbind", object$cf@prediction_weights(newdata = newdata, OOB = FALSE))
}
### method for for trtf::traforest
predict.mytrtf <- function(object, newdata = object$data, OOB = FALSE, type = "weights", ...)
{
    stopifnot(!OOB)
    predict(object$trtf, newdata = newdata, OOB = FALSE, type = "weights")
}
## method for randomForestSRC::rfsrc
predict.myrfsrc <- function(object, newdata = object$data, OOB = FALSE, type = "weights", ...)
{
    stopifnot(!OOB)
    nd <- newdata[,names(newdata) != "y"]   
    t(predict(object$rfsrc, newdata = nd, forest.wt = TRUE)$forest.wt)
}
### method for ranger::ranger
predict.myranger <- function(object, newdata = object$data, OOB = FALSE, 
                             type = "weights", ...) {
    stopifnot(!OOB)
    inbag <- object$rg$inbag.counts
    nnewdata <- predict(object$rg, data = newdata, type = "terminalNodes")$predictions
    ndata <- predict(object$rg, data = object$data, type = "terminalNodes")$predictions

    w <- 0L
    for (b in 1:length(inbag)) {
        fdata <- ndata[,b]
        fnewdata <- nnewdata[,b]
        ids <- unique(fdata)
        tw <- inbag[[b]]
        pw <- sapply(ids, function(i) tw * (fdata == i))
        ret <- pw[, match(fnewdata, ids), drop = FALSE]
        if (OOB) ret[,tw > 0] <- 0
        w <- w + ret
    }
    return(w)
}
### method for L1SRC::rfsrc
predict.myL1 <- function(object, newdata = object$data, OOB = FALSE, 
                         type = "weights", ...) {
    stopifnot(!OOB)
    inbag <- object$weights
    nnewdata <- predict(object$L1, newdata = newdata, membership = TRUE)$membership
    ndata <- predict(object$L1, newdata = object$data, membership = TRUE)$membership

    w <- 0L
    for (b in 1:length(inbag)) {
        fdata <- ndata[,b]
        fnewdata <- nnewdata[,b]
        ids <- unique(fdata)
        tw <- inbag[[b]]
        pw <- sapply(ids, function(i) tw * (fdata == i))
        ret <- pw[, match(fnewdata, ids), drop = FALSE]
        if (OOB) ret[,tw > 0] <- 0
        w <- w + ret
    }
    return(w)
}

### Transfromation survival forest with Bernstein basis and general score
### Bs(theta)
mytraforest <- function(learn, test, mtry = mtry, order = ord, log_first = logf) {
    
    ### the largest values of time are not used for Bernstein approximation
    sup <- c(sqrt(.Machine$double.eps), quantile(learn$y[,1], .95))
    m0 <- as.mlt(Coxph(y ~ 1, 
                       data = learn, 
                       order = order, 
                       log_first = log_first, 
                       sup = sup))

    ### traforest generates the subsampling weights later to be used
    ### by other methods
    rf <- traforest(m0, 
                    formula = y ~ ., 
                    data = learn,
                    ntree = ntree, 
                    trace = FALSE, 
                    mtry = mtry,
                    control = ctrl_partykit)

    ### nearest neighbor predictions
    cf <- predict(rf, newdata = test, type = "coef")
    llNN <- logLik(rf, newdata = test, coef = cf)
	
    list(logLik_NN = llNN, m0 = m0, object = rf, 
         mtry = mtry)

}

### Weibull distributional survival forests with general score
### W(theta)
mytraforest_Seibold <- function(trtf, learn, test) {
    mtry <- trtf$mtry
    m0 <- trtf$m0
    trtf <- trtf$object

    ### mytraforest generated the subsampling weights 
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"

    m0 <- as.mlt(Survreg(y ~ 1, data = learn))

    ### split wrt intercept and log(y)
    rf <- traforest(m0, 
                    formula = y ~ ., 
                    data = learn,
                    ntree = ntree, 
                    trace = FALSE, 
                    mtry = mtry,
                    control = ctrl_partykit, 
                    weights = initw)

    ### use hack to wire rf  weights into predict.traforest
    trtf$trtf <- rf
    class(trtf) <- c(class(trtf)[1], "mytrtf", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(logLik_NN = llNN, m0 = m0, object = rf)
}


### Weibull distributional survival forests with log-rank score
### W(alpha)
mytraforest_W <- function(trtf, learn, test) {                                                                                                                                     
    
    mtry <- trtf$mtry
    m0 <- trtf$m0
    trtf <- trtf$object

    ### mytraforest generated the subsampling weights 
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"
    
    m0 <- as.mlt(Survreg(y ~ 1, data = learn))

    ### split wrt to intercept only
    rf <- traforest(m0,
                    formula = y ~ ., 
                    data = learn, 
                    parm = "(Intercept)",
                    ntree = ntree, 
                    trace = FALSE, 
                    mtry = mtry,
                    control = ctrl_partykit, 
                    weights = initw)

    ### use hack to wire rf weights into predict.traforest
    trtf$trtf <- rf
    class(trtf) <- c(class(trtf)[1], "mytrtf", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(logLik_NN = llNN, m0 = m0, object = rf)

}

### Trnasformation survival forests with Bernstein basis and log-rank score
### Bs(alpha)
mytraforest_B <- function(trtf, learn, test, order = ord, log_first = logf) {
    mtry <- trtf$mtry
    m0 <- trtf$m0
    trtf <- trtf$object

    ### mytraforest generated the subsampling weights 
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"
    
    ### There is no intercept in Bernstein-basis, therefore to perform
    ### log-rank splitting constant variable is added to the data implicitly
    learn$alpha <- 1
    test$alpha <- 1
    
    ### the largest values of time are not used for Bernstein approximation
    sup <- c(sqrt(.Machine$double.eps), quantile(learn$y[, 1], .95))
    m0 <- as.mlt(Coxph(y ~ alpha, 
                       data = learn, 
                       order = order, 
                       log_first = log_first,
                       support = sup, 
                       fixed = c("alpha" = 0)))

    ### split wrt to intercept ("log-rank scores") only
    rf <- traforest(m0, 
                    formula = y | alpha ~ ., 
                    data = learn,
                    parm = "alpha",
                    mltargs = list(fixed = c("alpha" = 0)),
                    ntree = ntree, 
                    trace = FALSE, 
                    mtry = mtry,
                    control = ctrl_partykit, 
                    weights = initw)

    ### use hack to wire rf weights into predict.traforest
    trtf$trtf <- rf
    class(trtf) <- c(class(trtf)[1], "mytrtf", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(logLik_NN = llNN, m0 = m0, object = rf)
}


### cforests for survival data
mycforest <- function(trtf, learn, test) {

    mtry <- trtf$mtry
    m0 <- trtf$m0
    trtf <- trtf$object

    ctrl_party <- party:::cforest_unbiased(ntree = ntree, 
                                           maxdepth = nodedepth, 
                                           mtry = mtry, minbucket = minbucket)

    ### mytraforest generated the subsampling weights 
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"

    ### log-rank splitting
    rf <- party::cforest(y ~ ., data = learn, weights = initw, 
                         control = ctrl_party)

    ### use hack to wire rf weights into predict.traforest
    trtf$cf <- rf
    class(trtf) <- c(class(trtf)[1], "mycf", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)
    
    list(logLik_NN = llNN, object = rf)

}

### random forests SRC
myrfsrc <- function(trtf, learn, test) {

    mtry <- trtf$mtry
    m0 <- trtf$m0
    trtf <- trtf$object

    ### mytraforest generated the subsampling weights 
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"

    ### rfsrc needs Surv() in the formula, otherwise
    ### y is treated as multivariate response
    learn$time <- learn$y[,1]
    learn$event <- learn$y[,2]
    learn$y <- NULL
    ### use log-rank splitting
    rf <- randomForestSRC::rfsrc(Surv(time, event) ~ ., 
                                 data = learn, 
                                 samp = initw, 
                                 bootstrap = "by.user", 
                                 ntree = ntree, 
                                 nodedepth = nodedepth, 
                                 nodesize = minbucket, 
                                 mtry = mtry, 
                                 splitrule = "logrank")
    learn$y <- with(learn, Surv(time, event))

    ### use hack to wire rf weights into predict.traforest
    trtf$rfsrc <- rf
    class(trtf) <- c(class(trtf)[1], "myrfsrc", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(logLik_NN = llNN, object = rf)
}

### L1 splitting
myL1 <- function(trtf, learn, test) {

    mtry <- trtf$mtry
    m0 <- trtf$m0
    trtf <- trtf$object

    ### mytraforest generated the subsampling weights 
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"

    ### L1 needs Surv() in the formula, otherwise
    ### y is treated as multivariate response
    learn$time <- learn$y[,1]
    learn$event <- learn$y[,2]
    learn$y <- NULL
    ### use log-rank splitting
    rf <- L1SRC::rfsrc(Surv(time, event) ~ ., 
                       data = learn, 
                       samp = initw, 
                       bootstrap = "by.user", 
                       ntree = ntree, 
                       nodedepth = nodedepth, 
                       nodesize = minbucket, 
                       mtry = mtry, 
                       splitrule = "custom2") ### custom2 is L1 by Denis!
    learn$y <- with(learn, Surv(time, event))

    ### use hack to wire rf weights into predict.traforest
    trtf$L1 <- rf
    class(trtf) <- c(class(trtf)[1], "myL1", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)
    
    list(logLik_NN = llNN, object = rf)
}

### ranger
myranger <- function(trtf, learn, test) {

    mtry <- trtf$mtry
    m0 <- trtf$m0
    trtf <- trtf$object

    ### mytraforest generated the subsampling weights 
    initw <- do.call("cbind", trtf$weights)
    storage.mode(initw) <- "double"

    ### cannot use subsampling weights from mytraforest, cannot restrict
    ### maxdepths (use approximation via min.node.size)
    rf <- ranger(y ~ ., 
                 data = learn, 
                 num.trees = ntree, 
                 mtry = mtry, 
                 replace = FALSE, 
                 splitrule = "logrank",  
                 min.node.size = max(minbucket, nrow(learn) / 2^nodedepth), 
                 keep.inbag = TRUE)

    ### use hack to wire rf weights into predict.traforest
    trtf$rg <- rf
    class(trtf) <- c(class(trtf)[1], "myranger", class(trtf)[-1])

    ### nearest neighbor predictions
    cf <- predict(trtf, newdata = test, type = "coef")
    llNN <- logLik(trtf, newdata = test, coef = cf)

    list(logLik_NN = llNN, object = rf)

}

###simple example
if (FALSE) {
data("GBSG2", package = "TH.data")
GBSG2$y <- with(GBSG2, Surv(time, cens))
GBSG2$time <- GBSG2$cens <- NULL

ntree <- 25
test <- sample(1:NROW(GBSG2), 25)

tf <- mytraforest(GBSG2[-test,], GBSG2[test,], mtry = 3)
tfS <- mytraforest_Seibold(tf, GBSG2[-test,], GBSG2[test,])
tfB <- mytraforest_B(tf, GBSG2[-test,], GBSG2[test,])
tfW <- mytraforest_W(tf, GBSG2[-test,], GBSG2[test,])
cf <- mycforest(tf, GBSG2[-test,], GBSG2[test,])
library("randomForestSRC")
sf <- myrfsrc(tf, GBSG2[-test,], GBSG2[test,])
library("ranger")
rg <- myranger(tf, GBSG2[-test,], GBSG2[test,])
##not public; private patch from Denis Laroque
detach(package:randomForestSRC)
library("L1SRC")
L1 <- myL1(tf, GBSG2[-test,], GBSG2[test,])
detach(package:L1SRC)

tf[1]
tfS[1]
tfB[1]
tfW[1]
cf[1]
sf[1]
rg[1]
L1[1]

}
