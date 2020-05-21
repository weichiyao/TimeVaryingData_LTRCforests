
.is.formula <- function(x) {
    if (is.null(x)) return(FALSE)
    inherits(x, "formula")
}

.is.Surv <- function(x)
    inherits(x, "Surv")

.is.R <- function(x)
    inherits(x, "response")

.type_of_response <- function(y) {
    if (.is.Surv(y)) return("survival")
    if (.is.R(y)) {
        if (any(.exact(y))) {
            y <- y$exact[.exact(y)]
        } else {
            y <- y$cleft[.cleft(y)]
        }
    }
    if (storage.mode(y) == "double") return("double")
    if (is.integer(y)) return("integer")
    if (is.ordered(y)) return("ordered")
    if (is.factor(y)) return("unordered")
    return(NA)
}

.checkR <- function(x, weights = NULL) {
    if (!.is.R(x)) return(FALSE)
    if (is.null(weights)) weights <- 1
    if (all(.cleft(x) & !.cright(x) & weights > 0)) {
        warning("response contains left-censored observations only")
        return(FALSE)
    }
    if (all(!.cleft(x) & .cright(x) & weights > 0)) {
        warning("response contains right-censored observations only")
        return(FALSE)
    }
    return(TRUE)
}

.exact_subset <- function(exact, subset = NULL) {

    iex <- inex <- NULL
    if (any(exact))
        iex <- which(exact)
    if (any(!exact))
        inex <- which(!exact)
    if (is.null(subset))
        return(list(full_ex = iex, full_nex = inex))

    full_ex <- redu_ex <- full_nex <- redu_nex <- NULL
    if (length(iex) > 0) {
        full_ex <- iex[iex %in% subset]
        redu_ex <- (1L:length(iex))[iex %in% subset]
    }
   if (length(inex) > 0) {
        full_nex <- inex[inex %in% subset]
        redu_nex <- (1L:length(inex))[inex %in% subset]
    }
    return(list(full_ex = full_ex, redu_ex = redu_ex,
                full_nex = full_nex, redu_nex = redu_nex))
}

if (FALSE) {
exact <- c(TRUE, FALSE)[c(1, 2, 1, 1, 2, 2)]
.exact_subset(rep(TRUE, 3))
.exact_subset(rep(FALSE, 3))
.exact_subset(exact)
.exact_subset(exact, 1:2)
.exact_subset(exact, 3:6)
}
