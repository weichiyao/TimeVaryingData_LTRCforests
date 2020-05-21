\name{tram-methods}
\alias{as.mlt.tram}
\alias{model.frame.tram}
\alias{model.matrix.tram}
\alias{coef.tram}
\alias{coef.Lm}
\alias{coef.Survreg}
\alias{vcov.tram}
\alias{logLik.tram}
\alias{estfun.tram}
\alias{predict.tram}
\alias{plot.tram}
\title{
  Methods for Stratified Linear Transformation Models
}
\description{
  Methods for objects inheriting from class tram
}
\usage{
\method{as.mlt}{tram}(object)
\method{model.frame}{tram}(formula, ...)
\method{model.matrix}{tram}(object, data = object$data, with_baseline = FALSE, ...) 
\method{coef}{tram}(object, with_baseline = FALSE, ...) 
\method{coef}{Lm}(object, as.lm = FALSE, ...)
\method{coef}{Survreg}(object, as.survreg = FALSE, ...)
\method{vcov}{tram}(object, with_baseline = FALSE, ...) 
\method{logLik}{tram}(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
\method{estfun}{tram}(object, parm = coef(as.mlt(object), fixed = FALSE), ...)
\method{predict}{tram}(object, newdata = model.frame(object), 
        type = c("lp", "trafo", "distribution", "survivor", "density", 
                    "logdensity", "hazard", "loghazard", "cumhazard", "quantile"), 
           ...)
\method{plot}{tram}(x, newdata = model.frame(x), 
     which = c("QQ-PIT", "baseline only", "distribution"), 
     confidence = c("none", "interval", "band"), level = 0.95, 
     K = 50, cheat = K, col = "black", fill = "lightgrey", lwd = 1, ...)
}
\arguments{
  \item{object, formula, x}{a fitted stratified linear transformation model inheriting
                from class \code{tram}.}
  \item{data}{an optional data frame.}
  \item{with_baseline}{logical, if \code{TRUE} all model parameters
                       are returned, otherwise parameters describing the
                       baseline transformation are ignored.}
  \item{as.lm}{logical, return parameters in the \code{\link[stats]{lm}}
               parameterisation if \code{TRUE}.}
  \item{as.survreg}{logical, return parameters in the \code{\link[survival]{survreg}}
               parameterisation in \code{TRUE}.}
  \item{parm}{model parameters, including baseline parameters.}
  \item{newdata}{an optional data frame of new observations.}
  \item{type}{type of prediction, current options include
              linear predictors (\code{"lp"}, of \code{x} variables in the
              formula \code{y | s ~ x}), transformation functions
              (\code{"trafo"}) or distribution functions on the
              scale of the cdf (\code{"distribution"}),
              survivor function, density function, log-density
              function, hazard function, log-hazard function, cumulative
              hazard function or quantile function.}
  \item{which}{type of plot, either a QQ plot of the probability-integral
               transformed observations (\code{"QQ-PIT"}), of the
               baseline transformation of the whole distribution.}
  \item{confidence}{type of uncertainty assessment.}
  \item{level}{confidence level.}
  \item{K}{number of grid points in the response, see
           \code{\link[mlt]{plot.ctm}}.}
  \item{cheat}{reduced number of grid points for the computation
               of confidence bands, see \code{\link[mlt]{confband}}.}
  \item{col}{line color.}
  \item{fill}{fill color.}
  \item{lwd}{line width.}
  \item{\dots}{additional arguments to the underlying methods for class
               \code{mlt}, see \code{\link[mlt]{mlt-methods}}.}
}
\details{
  \code{coef} can be used to get (and set) model parameters, 
  \code{logLik} evaluates the log-likelihood (also for
  parameters other than the maximum likelihood estimate); 
  \code{vcov} returns the estimated variance-covariance matrix (possibly
  taking \code{cluster} into account) and 
  and \code{estfun} gives the score contribution by each observation.
  \code{predict} and \code{plot} can be used to inspect the model on
  different scales.
}
\seealso{
  \code{\link[mlt]{mlt-methods}},  \code{\link[mlt]{plot.ctm}}
}
\references{
  Torsten Hothorn, Lisa Moest, Peter Buehlmann (2018), Most Likely
  Transformations, \emph{Scandinavian Journal of Statistics}, \bold{45}(1),
  110--134, \doi{10.1111/sjos.12291}.
}
\examples{

    data("BostonHousing2", package = "mlbench")

    ### fit non-normal Box-Cox type linear model with two
    ### baseline functions (for houses near and off Charles River)
    BC_BH_2 <- BoxCox(cmedv | 0 + chas ~ crim + zn + indus + nox + 
                      rm + age + dis + rad + tax + ptratio + b + lstat,
                      data = BostonHousing2)
    logLik(BC_BH_2)

    ### classical likelihood inference
    summary(BC_BH_2)

    ### coefficients of the linear predictor
    coef(BC_BH_2)

    ### plot linear predictor (mean of _transformed_ response) 
    ### vs. observed values
    plot(predict(BC_BH_2, type = "lp"), BostonHousing2$cmedv)

    ### all coefficients
    coef(BC_BH_2, with_baseline = TRUE)

    ### compute predicted median along with 10% and 90% quantile for the first
    ### observations
    predict(BC_BH_2, newdata = BostonHousing2[1:3,], type = "quantile",
            prob = c(.1, .5, .9))

    ### plot the predicted density for these observations
    plot(BC_BH_2, newdata = BostonHousing2[1:3, -1],
         which = "distribution", type = "density", K = 1000)

    ### evaluate the two baseline transformations, with confidence intervals
    nd <- model.frame(BC_BH_2)[1:2, -1]
    nd$chas <- factor(c("0", "1"))
    library("colorspace")
    col <- diverge_hcl(2, h = c(246, 40), c = 96, l = c(65, 90))
    fill <- diverge_hcl(2, h = c(246, 40), c = 96, l = c(65, 90), alpha = .3)
    plot(BC_BH_2, which = "baseline only", newdata = nd, col = col,
         confidence = "interval", fill = fill, lwd = 2,
         xlab = "Median Value", ylab = expression(h[Y]))
    legend("bottomright", lty = 1, col = col, 
            title = "Near Charles River", legend = c("no", "yes"), bty = "n")

}
