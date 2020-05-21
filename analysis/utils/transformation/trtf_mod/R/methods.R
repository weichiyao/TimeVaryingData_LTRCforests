
.R2vec <- function(object) {
    if (!inherits(object, "response")) return(object)
    ex <- object$exact
    le <- object$cleft
    ri <- object$cright
    rex <- ex
    rle <- le
    rri <- ri
    rex[is.na(ex)] <- 0
    rle[is.na(le) | !is.finite(le)] <- 0
    rri[is.na(ri) | !is.finite(ri)] <- 0
    ### (-Inf, x] -> x and (x, Inf) -> x
    rex + (rle + ifelse(is.finite(ri) & is.finite(le), (rri - rle)/2, rri))
}


coef.trafotree <- function(object, ...)
    object$coef

logLik.trafotree <- function(object, newdata, weights = NULL, perm = NULL, ...) {

    cf <- coef(object)
    if (missing(newdata) && is.null(perm) && is.null(weights)) {
        ### <FIXME> this is different when nmax < Inf </FIXME>
        ret <- sum(object$logLik)
    } else {
        if (missing(newdata)) {
            newdata <- object$data
            if (is.null(weights))
                weights <- data_party(object)[["(weights)"]]
        } else {
            if (is.null(weights))
                weights <- rep(1, nrow(newdata))
        }
        tids <- nodeids(object, terminal = TRUE)
        nd <- factor(predict(object, newdata = newdata, type = "node", perm = perm, ...), 
                     levels = tids, labels = tids)
        ### set up unfitted model with newdata
        mltargs <- object$mltargs
        mltargs$data <- newdata
        mltargs$dofit <- FALSE
        mltmod <- do.call("mlt", mltargs)
        ndcf <- cf[unclass(nd),,drop = FALSE]
        if (NROW(ndcf) == 1L) ndcf <- ndcf[1,,drop = TRUE]
        ret <- logLik(mltmod, parm = ndcf, w = weights)
    }
    attr(ret, "df") <- length(cf)
    class(ret) <- "logLik"
    ret
}

logLik.traforest <- function(object, newdata, weights = NULL, OOB = FALSE, coef = NULL,  ...) {

    if (is.null(coef)) {
        if (missing(newdata)) {
            cf <- predict(object, OOB = OOB, type = "coef")
        } else {
            cf <- predict(object, newdata = newdata, type = "coef")
        }
    } else {
        cf <- coef
    }

    if (missing(newdata)) {
        newdata <- object$data
        if (is.null(weights))
            weights <- object$fitted[["(weights)"]]
    } 
    ### set up unfitted model with newdata
    mltargs <- object$mltargs
    mltargs$data <- newdata
    mltargs$dofit <- FALSE
    mltmod <- do.call("mlt", mltargs)

    if (is.null(weights)) weights <- rep(1, nrow(newdata))
    ndcf <- do.call("rbind", cf)
    if (NROW(ndcf) == 1L) ndcf <- ndcf[1,,drop = TRUE]
    ret <- logLik(mltmod, parm = ndcf, w = weights)
    attr(ret, "df") <- NA
    class(ret) <- "logLik"
    ret
}

predict.trafotree <- function(object, newdata, K = 20, q = NULL,
    type = c("node", "coef", "trafo", "distribution", "survivor", "density",
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"),
    perm = NULL, ...) {

    type <- match.arg(type)
    tmp <- object
    class(tmp) <- class(tmp)[-1L]
    if (missing(newdata)) {
        nf <- predict(tmp, type = "node", perm = perm)
    } else {
        nf <- predict(tmp, newdata = newdata, type = "node", perm = perm)
    }
    if (type == "node") return(nf)
    tids <- nodeids(object, terminal = TRUE)
    nf <- factor(nf, levels = tids, labels = tids)
    if (type == "coef") return(object$coef[nf,,drop = FALSE])

    mod <- object$model
    if (is.null(q))
        q <- mkgrid(mod, n = K)[[mod$response]]

    if (missing(newdata)) newdata <- data_party(object)

    ### <FIXME> need .R2vec??? </FIXME>
    pr <- .R2vec(predict(mod, newdata = newdata, q = q, type = type, ...))
    if (!is.matrix(pr))
        pr <- matrix(pr, nrow = NROW(pr), ncol = NROW(newdata))
    for (nd in levels(nf)) {
        i <- nf == nd
        coef(mod) <- object$coef[nd,]
        pr[,i] <- .R2vec(predict(mod, newdata = newdata[i,,drop = FALSE], q = q,
                                 type = type, ...))
    } 
    # ret<- NULL
    # ret$q = q
    # ret$ans = pr
    return(pr)
}

.thetastart <- function(object, weights, i, updatestart, cf) {

    if (!updatestart) return(coef(object, fixed = FALSE))
    ### not sure if better starting values are worth the time,
    ### so coef(object) is still the default
    m <- object$model$model
    w <- weights[,i]
    if (i <= 25) {
        theta <- coef(object) 
        if (names(m) == "bresponse" && 
            inherits(m$bresponse, "Bernstein_basis")) {
            ### Bernstein: theta_k = f(k / n)
            y <- object$response
            su <- variables::support(attr(m$bresponse, "variables"))[[1]]
            grid <- seq(from = su[1], to = su[2], length.out = length(coef(object)))
            prob <- attr(y, "prob")(w)(grid)
            prob <- pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, prob))
            theta <- object$model$todistr$q(prob)
        }
    } else {
        imin <- which.min(cs <- colSums((weights[, 1:(i - 1), drop = FALSE] - w)^2))
        theta <- cf[[imin]]
    }
    return(theta)
}

predict.traforest <- function(object,  newdata, mnewdata = data.frame(1), K = 20, q = NULL,
    type = c("weights", "node", "coef", "trafo", "distribution", "survivor", "density",
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"),
    OOB = FALSE, simplify = FALSE, trace = FALSE, updatestart = FALSE, 
    applyfun = NULL, cores = NULL, ...) {

    type <- match.arg(type)

    if (!missing(newdata) && !type == "node" && (!is.null(applyfun) || !is.null(cores))) {
        call <- match.call()
        if (is.null(applyfun)) {
            applyfun <- if(is.null(cores)) {
                lapply  
            } else {
                function(X, FUN, ...)
                    parallel::mclapply(X, FUN, ..., mc.cores = cores)
            }
        }
        i <- 1:nrow(newdata)
        idx <- cut(i, breaks = (0:cores/cores) * nrow(newdata))
        idx <- tapply(i, idx, function(x) x, simplify = FALSE)
        call$applyfun <- call$cores <- NULL
        ret <- applyfun(idx, function(i) {
            predict(object = object, newdata = newdata[i,,drop = FALSE],
                    mnewdata = mnewdata, K = K, q = q, type = type, 
                    OOB = OOB, simplify = simplify, trace = trace, 
                    updatestart = updatestart, applyfun = NULL, cores = NULL, ...)
        })
        type <- match.arg(type)
        names(ret) <- NULL
        if (type == "weights") {
            ret <- do.call("cbind", ret)
        } else {
            ret <- do.call("c", ret)
        }
        return(ret)
    }

    tmp <- object
    class(tmp) <- class(tmp)[-1L]

    ptype <- type
    if (!(ptype %in% c("weights", "node"))) ptype <- "weights"
    if (missing(newdata)) {
        ret <- predict(tmp, OOB = OOB, type = ptype,
                       simplify = TRUE)
    } else {
        ret <- predict(tmp, newdata = newdata, type = ptype,
                       simplify = TRUE)
    }
    if (type %in% c("weights", "node")) return(ret)

    mod <- object$model
    if (is.null(q))
        q <- mkgrid(mod, n = K)[[mod$response]]
    mltmod <- object$mltobj
    mod <- mltmod$object
    thetastart <- coef(mod, fixed = FALSE)

    ans <- vector(mode = "list", length = ncol(ret))
    names(ans) <- colnames(ret)
    cf <- vector(mode = "list", length = ncol(ret))
    names(cf) <- colnames(ret)

    converged <- logical(ncol(ret))
    if (trace) pb <- txtProgressBar(style = 3)
    for (i in 1:ncol(ret)) {
        if (trace) setTxtProgressBar(pb, i / ncol(ret))
        w <- ret[,i]
        thetastart <- .thetastart(mltmod$object, ret, i, updatestart, cf)
        ### try hard to estimate parameters; if may happen that parameters for 
        ### a specific obs are not identified (out of range)
        umod <- try(object$trafo(subset = which(w > 0), 
                                 weights = w, info = list(coef = thetastart),
                                 estfun = FALSE), silent = TRUE)
        converged[i] <- umod$converged 
        if (inherits(umod, "try-error")) {
            cf[[i]] <- NA
            ans[[i]] <- NA
        } else {
            cf[[i]] <- umod$coef
            if (type != "coef") {
                coef(mod) <- umod$coef
                ans[[i]] <- predict(mod, q = q, newdata = mnewdata, type = type, ...)
            }
        }
    } 
    if (!all(converged)) 
        warning("Parameter estimation did not converge for observations ",
                paste(which(!converged), sep = ", "))
    if (trace) close(pb)
    if (type == "coef") return(cf)
    # if (type == "survivor"){
    #   ret<-NULL
    #   ret$ans <- ans
    #   ret$q <- q
    # }
    return(ans)
}

simulate.traforest <- function(object, nsim = 1, seed = NULL, newdata, 
                               OOB = FALSE, coef = NULL, ...) {

    if (is.null(coef)) {
        if (missing(newdata)) {
            cf <- predict(object, type = "coef", OOB = OOB)
            newdata <- object$data
        } else {
            cf <- predict(object, newdata = newdata, type = "coef")
        }
    } else {
        cf <- coef
        newdata <- object$data
    }
    if (is.list(cf)) cf <- do.call("rbind", cf)
    if (nrow(cf) != nrow(newdata)) stop("coef and newdata don't match")

    mod <- object$model
    ret <- vector(mode = "list", length = nrow(cf))
    for (i in 1:nrow(cf)) {
        coef(mod) <- cf[i,]
        ret[[i]] <- simulate(mod, nsim = nsim, seed = seed, 
                             newdata = newdata[i,,drop = FALSE], bysim = FALSE, ...)[[1]]
    }
    ans <- vector(mode = "list", length = nsim)
    if (any(sapply(ret, function(x) inherits(x, "response")))) {
        ret <- lapply(ret, function(x) {
            if (inherits(x, "response")) return(x)
            R(x)
        })
        for (j in 1:nsim) {
            for (i in 1:nrow(cf))
                ans[[j]] <- rbind(ans[[j]], ret[[i]][j,])
        }
    } else {
        for (j in 1:nsim) {
            for (i in 1:nrow(cf))
                ans[[j]] <- rbind(ans[[j]], ret[[i]][j])
        }
    }
    ans
}

simulate.trafotree <- simulate.traforest

gettree.traforest <- function(object, tree = 1L, ...) {

    ft <- object$fitted
    ft[["(weights)"]] <- weights <- object$weights[[tree]]
    ret <- party(object$nodes[[tree]], data = object$data, fitted = ft)
    ret$terms <- object$terms
    class(ret) <- c("constparty", class(ret))
    ret$model <- object$model
    ret$mltobj <- object$mltobj
    ret$mltargs <- object$mltargs
    ret$trafo <- object$trafo
    nd <- predict(ret, newdata = object$data, type = "node")
    ret$models <- tapply(1:length(nd), factor(nd), function(i)
        ret$trafo(i, weights = weights, estfun = FALSE)) ### note: trafo needs weights
    ret$coef <- do.call("rbind", lapply(ret$models, function(x) x$coef))
    ret$logLik <- sapply(ret$models, function(x) -x$objfun)
    class(ret) <- c("trafotree", class(ret))
    ret
}
