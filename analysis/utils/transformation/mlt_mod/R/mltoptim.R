
mltoptim <- function(auglag = list(maxtry = 5, kkt2.check = FALSE), 
                     spg = list(maxit = 10000, quiet = TRUE, checkGrad = FALSE),
                     trace = FALSE) 
{
    ret <- list()
    if (!is.null(auglag))
        ret$auglag <- function(theta, f, g, ui = NULL, ci = NULL) {
            control <- auglag
            maxtry <- control$maxtry
            control$maxtry <- NULL
            control$trace <- trace
            atheta <- theta
            if (!is.null(ui)) {
                for (tr in 1:maxtry) {
                    ret <- try(alabama::auglag(par = atheta, fn = f, gr = g, hin = function(par) ui %*% par - ci, hin.jac = function(par) ui,
                                               control.outer = control)[c("par", "convergence", "value")])
                    atheta <- runif(length(atheta))
                    if (inherits(ret, "try-error")) next()
                    if (ret$convergence == 0) break()
                }
            } else { 
                control <- spg
                control$trace <- trace
                quiet <- control$quiet
                control$quiet <- NULL
                ret <- try(BBoptim(par = theta, fn = f, gr = g, control = control, quiet = quiet))
            }
            if (inherits(ret, "try-error"))
                ret <- list(par = theta, convergence = 1)
            return(ret)
        }
    if (!is.null(spg)) 
        ret$spg <- function(theta, f, g, ui = NULL, ci = NULL) {
            control <- spg
            control$trace <- trace
            quiet <- control$quiet
            control$quiet <- NULL
            if (!is.null(ui)) {
                ret <- try(BB::BBoptim(par = theta, fn = f, gr = g, project = "projectLinear",
                                       projectArgs = list(A = ui, b = ci, meq = 0), control = control,
                                       quiet = quiet))
            } else { 
                ret <- try(BBoptim(par = theta, fn = f, gr = g, control = control, quiet = quiet))
            }
            if (inherits(ret, "try-error"))
                ret <- list(par = theta, convergence = 1)
            return(ret)
        }
    return(ret)
}
