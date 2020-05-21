
library("trtf")
library("randomForest")
library("parallel")
RNGkind("L'Ecuyer-CMRG")
set.seed(29)
library("sandwich")

trafotree <- function(...) {
    ret <- try(trtf::trafotree(...))
    if (inherits(ret, "try-error")) return(NULL)
    return(ret)
}

traforest <- function(...) {
    ret <- try(trtf::traforest(...))
    if (inherits(ret, "try-error")) return(NULL)
    return(ret)
}


Ctree <- function(tt, root_trafo_only = FALSE, ...) {

    if (root_trafo_only) {
        ### NOTE: Setting up ytrafo that way ignores weights and thus we can't 
        ### use it for cforest !!!
        sc <- estfun(tt$trafo(model = TRUE)$object)
        ytrafo <- list(y = function(...) sc)
    } else {
        ytrafo <- NULL
    }
    ct <- ctree(..., ytrafo = ytrafo)
    ct$model <- tt$model
    ct$mltobj <- tt$trafo(model = TRUE, estfun = FALSE)
    ct$trafo <- tt$trafo
    ct$mltargs <- tt$mltargs

    nd <- predict(ct, type = "node")
    ct$models <- tapply(1:length(nd), factor(nd), function(i) 
        ct$trafo(i, estfun = FALSE)) ### note: trafo is (potentially) weighted
    ct$coef <- do.call("rbind", lapply(ct$models, function(x) x$coef))
    ct$logLik <- sapply(ct$models, function(x) x$objfun)

    class(ct) <- c("trafotree", class(ct))
    ct
}

brf <- function(...)
    ### NOTE: nodesize is minsplit NOT minbucket (in rpart terminology)
    randomForest(..., keep.inbag = TRUE, replace = FALSE, nodesize = 25)

myrf <- function(rf, trtf) {
    trtf$rf <- rf
    class(trtf) <- c(class(trtf)[1], "myrf", class(trtf)[-1])
    trtf
}

mycf <- function(cf, trtf) {
    trtf$cf <- cf
    class(trtf) <- c(class(trtf)[1], "mycf", class(trtf)[-1])
    trtf
}

predict.mycf <- function(object, newdata = object$data, OOB = FALSE, type = "weights", ...)
    predict(object$cf, newdata = newdata, OOB = FALSE, type = "weights")

predict.myrf <- function(object, newdata = object$data, OOB = FALSE, type = "weights", ...) {

    type <- match.arg(type)

    lnodes <- attr(predict(object$rf, newdata = object$data, nodes = TRUE), "nodes")

    fdata <- as.list(as.data.frame(lnodes))
    if (OOB) {
        fnewdata <- list()
    } else {
        nnodes <- attr(predict(object$rf, newdata = newdata, nodes = TRUE), "nodes")
        fnewdata <- as.list(as.data.frame(nnodes))
    }

    ### extract weights
    rw <- as.list(as.data.frame(object$rf$inbag))

    w <- partykit:::.rfweights(fdata, fnewdata, rw)

    return(w)
}
    
simfun <- function(..., nsim = 100, methods = c("ct", "tt", "ttExSplit", "ttM", "ttSD", "ttBern", 
    "ttBernExSplit", "tfbagg", "tfrf", "tfrfExSplit", "tfM", "tfSD", "tfbaggBern", "tfrfBern", 
    "tfrfBernExSplit", "cfbagg", "cfrf", "bagg", "rf", "rfBern"), prob = c(.1, .5, .9), bounds = c(-Inf, Inf)) {

    methods <- match.arg(methods, several.ok = TRUE)

    ll <- matrix(NA, nrow = nsim, ncol = length(methods) + 1)
    colnames(ll) <- c("gt", methods)
    time <- matrix(NA, nrow = nsim, ncol = length(methods) + 1)
    colnames(time) <- colnames(ll)
    covr <- matrix(NA, nrow = nrow(ll), ncol = ncol(ll))
    colnames(covr) <- colnames(ll)
    clen <- matrix(NA, nrow = nrow(ll), ncol = ncol(ll))
    colnames(clen) <- colnames(ll)
    c_covr <- matrix(NA, nrow = nrow(ll), ncol = ncol(ll))
    colnames(c_covr) <- colnames(ll)
    c_clen <- matrix(NA, nrow = nrow(ll), ncol = ncol(ll))
    colnames(c_clen) <- colnames(ll)
    L1 <- matrix(NA, nrow = nrow(ll), ncol = ncol(ll))
    colnames(L1) <- colnames(ll)
    checkRisk <- array(NA, dim = c(nrow(ll), ncol(ll), length(prob)),
                       dimnames = list(NULL, colnames(ll), as.character(prob)))

for (i in 1:nrow(ll)) {

    print(i)
    l1 <- dgp(...)
    args <- list(...)
    args$n <- args$ntest
    t1 <- do.call("dgp", args)

    fm <- y ~ .
    mtrybagg <- max(1, ncol(l1) - 1)
    mtryrf <- max(1, floor((ncol(l1) - 1) / 3))

    b <- as.basis(~ y, data = l1, ui = matrix(0:1, nrow = 1), ci = 0)
    m <- ctm(b, todistr = "Normal")

    qu <- quantile(l1$y, prob = c(.05, .95))
    add <- c(min(l1$y) - qu[1], max(l1$y) - qu[2])
    if (is.finite(bounds[1])) add[1] <- 0
    yvar <- numeric_var("y", support = qu, bounds = bounds, add = add)
    mBern <- ctm(Bernstein_basis(yvar, ui = "increasing", order = 5), todistr = "Normal")

    models <- vector(mode = "list", length = length(methods))
    names(models) <- methods

    ctrlExSplit <- ctrl
    ctrlExSplit$splitflavour <- "exhaustive"
    ctrlExSplit$restart <- FALSE

    tm <- function(...) system.time(...)[["user.self"]]

    if ("tt" %in% methods)
        time[i, "tt"] <- tm(models[["tt"]] <- trafotree(m, parm = 1:2, formula = fm, data = l1))

    if ("ct" %in% methods) {
       stopifnot(!is.null(models[["tt"]]))
       time[i, "ct"] <- tm(models[["ct"]] <- Ctree(formula = fm, data = l1, tt = models[["tt"]]))
    }

    if ("ttExSplit" %in% methods) {
        ttctrlExSplit <- ctree_control()
        ttctrlExSplit$splitflavour <- "exhaustive"
        ttctrlExSplit$restart <- FALSE
        time[i, "ttExSplit"] <- tm(models[["ttExSplit"]] <- trafotree(m, parm = 1:2, formula = fm, 
            data = l1, control = ttctrlExSplit))
    }

    if ("ttM" %in% methods)
        time[i, "ttM"] <- tm(models[["ttM"]] <- trafotree(m, parm = 1, formula = fm, data = l1))

    if ("ttSD" %in% methods)
        time[i, "ttSD"] <- tm(models[["ttSD"]] <- trafotree(m, parm = 2, formula = fm, data = l1))

    if ("tfbagg" %in% methods)
        time[i, "tfbagg"] <- tm(models[["tfbagg"]] <- traforest(m, parm = 1:2, formula = fm, data = l1, mtry = mtrybagg,
                        ntree = ntree, control = ctrl))

    if ("tfrf" %in% methods)
        time[i, "tfrf"] <- tm(models[["tfrf"]] <- traforest(m, parm = 1:2, formula = fm, data = l1, mtry = mtryrf,
                        ntree = ntree, control = ctrl))

    if ("tfrfExSplit" %in% methods)
        time[i, "tfrfExSplit"] <- tm(models[["tfrfExSplit"]] <- traforest(m, parm = 1:2, formula = fm, data = l1, mtry = mtryrf,
                        ntree = ntree, control = ctrlExSplit))

    if ("cfbagg" %in% methods) {
        stopifnot(!is.null(models[["tfbagg"]]))
        time[i, "cfbagg"] <- tm(models[["cfbagg"]] <- mycf(cforest(formula = fm, data = l1, mtry = mtrybagg,
                        ntree = ntree, control = ctrl), trtf = models[["tfbagg"]]))
    }

    if ("cfrf" %in% methods) {
        stopifnot(!is.null(models[["tfrf"]]))
        time[i, "cfrf"] <- tm(models[["cfrf"]] <- mycf(cforest(formula = fm, data = l1, mtry = mtryrf,
                        ntree = ntree, control = ctrl), trtf = models[["tfrf"]]))
    }

    if ("tfbaggM" %in% methods)
        time[i, "tfbaggM"] <- tm(models[["tfbaggM"]] <- traforest(m, parm = 1, formula = fm, data = l1, mtry = mtrybagg,
                                     ntree = ntree, control = ctrl))

    if ("tfbaggSD" %in% methods)
        time[i, "tfbaggSD"] <- tm(models[["tfbaggSD"]] <- traforest(m, parm = 2, formula = fm, data = l1, mtry = mtrybagg,
                                      ntree = ntree, control = ctrl))

    if ("ttBern" %in% methods)
        time[i, "ttBern"] <- tm(models[["ttBern"]] <- trafotree(mBern, formula = fm, data = l1))

    if ("ttBernExSplit" %in% methods) {
        ttctrlExSplit <- ctree_control()
        ttctrlExSplit$splitflavour <- "exhaustive"
        ttctrlExSplit$restart <- FALSE
        time[i, "ttBernExSplit"] <- tm(models[["ttBernExSplit"]] <- trafotree(mBern, formula = fm, 
            data = l1, control = ttctrlExSplit))
    }

    if ("tfbaggBern" %in% methods)
        time[i, "tfbaggBern"] <- tm(models[["tfbaggBern"]] <- traforest(mBern, formula = fm, data = l1, mtry = mtrybagg,
                            ntree = ntree, control = ctrl))

    if ("tfrfBern" %in% methods)
        time[i, "tfrfBern"] <- tm(models[["tfrfBern"]] <- traforest(mBern, formula = fm, data = l1, mtry = mtryrf,
                                          ntree = ntree, control = ctrl))

    if ("tfrfBernExSplit" %in% methods)
        time[i, "tfrfBernExSplit"] <- tm(models[["tfrfBernExSplit"]] <- traforest(mBern, formula = fm, data = l1, mtry = mtryrf,
                                          ntree = ntree, control = ctrlExSplit))

    if ("bagg" %in% methods) {
        stopifnot(!is.null(models[["tfbagg"]]))
        time[i, "bagg"] <- tm(models[["bagg"]] <- myrf(brf(fm, data = l1, ntree = ntree, mtry = mtrybagg), trtf = models[["tfbagg"]]))
    }

    if ("rf" %in% methods) {
        stopifnot(!is.null(models[["tfrf"]]))
        time[i, "rf"] <- tm(models[["rf"]] <- myrf(brf(fm, data = l1, ntree = ntree, mtry = mtryrf), trtf = models[["tfrf"]]))
    }

    if ("rfBern" %in% methods && !is.null(models[["tfrfBern"]])) {
        stopifnot(!is.null(models[["tfrfBern"]]))
        time[i, "rfBern"] <- tm(models[["rfBern"]] <- myrf(brf(fm, data = l1, ntree = ntree, mtry = mtryrf), trtf = models[["tfrfBern"]]))
    }

    ll[i, "gt"] <- loglik(t1, ...)

    for (m in methods) {
        if (is.null(models[[m]])) next()
        ll[i,m] <- logLik(models[[m]], newdata = t1)
    }

    tmp <- t1
    tmp$y <- NULL

    quant <- vector(mode = "list", length = length(methods) + 1L)
    names(quant) <- c("gt", methods)
    quant[["gt"]] <- qu(t1, prob = prob, ...)

    for (m in methods) {
        if (is.null(models[[m]])) next()
        qtmp <- predict(models[[m]], newdata = tmp, type = "quantile", prob = prob, K = 500)
        if (is.list(qtmp)) {
            qtmp <- lapply(qtmp, trtf:::.R2vec)
            qtmp <- do.call("cbind", qtmp)
        }
        if (is.matrix(qtmp)) qtmp <- t(qtmp)
        quant[[m]] <- qtmp
    }

    for (m in c("gt", methods)) {
        if (m != "gt") {
            if (is.null(models[[m]])) next()
        }
        a <- try(u <- t1$y - quant[[m]])
        a <- try(um <- abs(u) * (u < 0) * matrix(1 - prob, ncol = length(prob), nrow = NROW(u), byrow = TRUE))
        a <- try(up <- abs(u) * (u >= 0) * matrix(prob, ncol = length(prob), nrow = NROW(u), byrow = TRUE))
        a <- try(checkRisk[i,m,] <- colSums(um + up))
    }
}

return(list(parm = list(...), ll = ll, time = time, checkRisk = checkRisk))
}
