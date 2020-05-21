
ctm <- function(response, interacting = NULL, shifting = NULL, data = NULL,
                todistr = c("Normal", "Logistic", "MinExtrVal", "MaxExtrVal"),
                sumconstr = inherits(interacting, c("formula", "formula_basis")), ...) {

    ### mkgrid() will not work if data is missing
    if (.is.formula(response)) 
        response <- as.basis(response, data = data)
    if (.is.formula(interacting)) 
        interacting <- as.basis(interacting, data = data)
    if (.is.formula(shifting)) 
        shifting <- as.basis(shifting, data = data, remove_intercept = TRUE, ...)

    if (is.character(todistr))
        todistr <- .distr(todistr)

    bases <- list(response = response, interacting = interacting, shifting = shifting)

    if (!is.null(interacting))
        interacting <- b(iresponse = response, iinteracting = interacting, 
                         sumconstr = sumconstr)

    if (is.null(interacting) && is.null(shifting)) {
        mod <- c(bresponse = response)
    } else if (!is.null(interacting) && is.null(shifting)) {
        mod <- c(binteracting = interacting)
    } else if (is.null(interacting) && !is.null(shifting)) {
        mod <- c(bresponse = response, bshifting = shifting)
    } else {
        mod <- c(binteracting = interacting, bshifting = shifting)
    }
    ret <- list(model = mod, response = variable.names(response), 
                todistr = todistr, bases = bases)
    class(ret) <- "ctm"
    nd <- lapply(mkgrid(ret, n = 1), function(x) x[1]) ### integer may have more values
    nd <- do.call("expand.grid", nd)
    X <- model.matrix(ret, data = nd)
    cf <- numeric(NCOL(X))
    names(cf) <- colnames(X)
    cf[] <- NA
    ret$coef <- cf
    return(ret)
}

model.matrix.ctm <- function(object, data, ...)
    return(model.matrix(object$model, data = data, ...))

variable.names.ctm <- function(object, 
    which = c("all", "response", "interacting", "shifting"), 
    ...) {

    which <- match.arg(which)
    m <- object$bases
    if (which == "all")
        return(variable.names(object$model))
    if (which == "response")
        return(variable.names(m$response))
    if (which == "interacting") {
        if (!is.null(m$interacting))
            return(variable.names(m$interacting))
        return(NULL)
    }
    if (which == "shifting") {
        if (!is.null(m$shifting))
            return(variable.names(m$shifting))
        return(NULL)
    }
}

coef.ctm <- function(object, ...)
    object$coef

"coef<-.ctm" <- function(object, value) {
    cf <- coef(object)
    stopifnot(length(cf) == length(value))
    if (!is.null(names(value))) {
        stopifnot(all(names(value) %in% names(cf)))
    } else {
        stopifnot(length(value) == length(cf))
        names(value) <- names(cf)
    }
    object$coef[names(value)] <- value
    object
}
