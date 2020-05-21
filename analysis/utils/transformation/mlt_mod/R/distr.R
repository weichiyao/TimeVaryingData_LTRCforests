
.Normal <- function()
    list(p = pnorm, d = dnorm, q = qnorm, 
         ### see also MiscTools::ddnorm
         dd = function(x) -dnorm(x = x) * x,
         ddd = function(x) dnorm(x = x) * (x^2 - 1), 
         dd2d = function(x) -x,
         name = "normal")

.Logistic <- function()
    list(p = plogis, d = dlogis, q = qlogis,
         dd = function(x) {
             ex <- exp(x)
             (ex - ex^2) / (1 + ex)^3
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 4 * ex^2 + ex^3) / (1 + ex)^4
         },
         dd2d = function(x) {
             ex <- exp(x)
             (1 - ex) / (1 + ex)
         },
         name = "logistic")

.MinExtrVal <- function()
    list(p = function(x) 1 - exp(-exp(x)),
         q = function(p) log(-log1p(- p)),
         d = function(x, log = FALSE) {
             ret <- x - exp(x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x) {
             ex <- exp(x)
             (ex - ex^2) / exp(ex)
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 3*ex^2 + ex^3) / exp(ex)
         },
         dd2d = function(x)
             1 - exp(x),
         name = "minimum extreme value")

.MaxExtrVal <- function()
    list(p = function(x) exp(-exp(-x)),
         q = function(p) -log(-log(p)),
         d = function(x, log = FALSE) {
             ret <- - x - exp(-x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x) {
             ex <- exp(-x)
             exp(-ex - x) * (ex - 1)
         },
         ddd = function(x) {
             ex <- exp(-x)
             exp(-x - ex) * (ex - 1)^2 - exp(-ex - 2 * x)
         },
         dd2d = function(x)
             exp(-x) - 1,
         name = "maximum extreme value")


.distr <- function(which = c("Normal", "Logistic", 
                             "MinExtrVal", "MaxExtrVal")) {
    which <- match.arg(which)
    do.call(paste(".", which, sep = ""), list())
}

