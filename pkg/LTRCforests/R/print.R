#' Print Summary Output of a ltrccif object or a ltrcrrf object
#'
#' Print summary output after a LTRCCIF or a LTRCRRF model is built.
#' This is the default print method for objects in the class of \code{\link{ltrccif}} or 
#' \code{\link{ltrcrrf}}.
#'
#' @param x an object of class \code{\link{ltrccif}} or \code{\link{ltrcrrf}}.
#' @return A printout object containing the following components:
#'         \item{Number of (pseudo-subject) observations}{number of left-truncated 
#'         right-censored pseudo-subject observations based on the Andersen-Gill reformulation.}
#'         \item{Number of subjects}{number of independent subject observations.}
#'         \item{Number of deaths}{number of times that an event occurs in the whole dataset.}
#'         \item{Number of trees}{the value set for argument \code{ntree}, 
#'         see \code{\link{ltrccif}} and \code{\link{ltrcrrf}}.}
#'         \item{minsplit}{the value set for argument \code{minsplit} that controls 
#'         the growth of individual trees; see \code{\link[partykit]{ctree_control}}.}
#'         \item{minbucket}{the value set for argument \code{minbucket} 
#'         that controls the growth of individual trees; see \code{\link[partykit]{ctree_control}}.}
#'         \item{minprob}{the value set for argument \code{minprob} 
#'         that controls the growth of individual trees; see \code{\link[partykit]{ctree_control}}.}
#'         \item{maxdepth}{the value set for argument \code{maxdepth} 
#'         that controls the maximum depth of individual trees; see \code{\link[partykit]{ctree_control}}.}
#'         \item{No. of variables tried at each split}{number of input variables 
#'         randomly sampled as candidates at each node for random forest algorithms, 
#'         which is either set as an argument \code{mtry} in \code{\link{ltrccif}} and \code{\link{ltrcrrf}}, 
#'         or tuned by \code{\link{tune.ltrccif}} or \code{\link{tune.ltrcrrf}}, respectively.}
#'         \item{Total no. of variables}{the number of features provided in \code{data}.}
#'         \item{Bootstrap type to grow trees}{the values set for augument \code{bootstrap}, 
#'         see \code{\link{ltrccif}} and \code{\link{ltrcrrf}}.}
#'         \item{Resampling used to grow trees}{the value set for argument \code{samptype}, 
#'         see \code{\link{ltrccif}} and \code{\link{ltrcrrf}}.}
#'         \item{Resampling rate used to grow trees}{the values set for argument \code{sampfrac}, 
#'         see \code{\link{ltrccif}} and \code{\link{ltrcrrf}}.}
#'         \item{Analysis}{LTRCCIF for a \code{\link{ltrccif}} object or LTRCRRF for \code{\link{ltrcrrf}}.}
#'         \item{Family}{the family used in the analysis, \code{surv}.}
#'         \item{Splitting rule}{the splitting rule that is implemented, 
#'         conditional inference framework for a \code{\link{ltrccif}} object or 
#'         Poisson splitting for \code{\link{ltrcrrf}}.}
#'         \item{Number of random split points}{the values set for argument \code{nsplit} in \code{\link{ltrcrrf}}.}
#' @importFrom survival Surv
#' @examples
#'
#' library(survival)
#' Formula = Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' # Built a LTRCCIF forest on the time-varying data by specifying id, with mtry specified:
#' LTRCCIFobj = ltrccif(formula = Formula, data = pbcsample, id = ID, mtry = 3, ntree = 50L)
#' print(LTRCCIFobj)
#'
#' # Built a LTRCCIF forest on the time-invariant data, with resampling, with mtry specified:
#' LTRCCIFobj = ltrccif(formula = Formula, data = pbcsample, samptype = "swr",
#'                      mtry = 3, ntree = 50L)
#' print(LTRCCIFobj)
#' @aliases print.ltrccif, print.ltrcrrf
#' @seealso \code{\link{ltrccif}}, \code{\link{ltrcrrf}}
#' @export

print <- function(x){
  UseMethod("print", x)
}

#' @export
print.ltrccif <- function(x) {
  ## check that the object is interpretable
  # if (any(inherits(x, c("ltrccif", "grow"), TRUE) ) == 0 ) {
  #   stop("This function only works for objects (of subclass) of class `(ltrccif, grow).")
  # }
  
  ## x will be processed if it's multivariate - therefore save some values from it
  familyPretty <- "LTRCCIF"
  familyOrg <- "surv"
  
  ## error rates
  err.rate <- NULL
  if (x$info$samptype == "swr" ){
    samptype = "sampling with replacement"
  } else {
    samptype = "sampling without replacement"
  }
  
  #################################################################################
  ##
  ## grow mode
  ##
  #################################################################################
  
  cat(" Number of (pseudo-subject) observations: ", nrow(x$fitted),                 "\n", sep="")
  
  cat("                      Number of subjects: ", length(unique(x$data$id)),  "\n", sep="")
  cat("                        Number of deaths: ", sum(as.matrix(x$fitted)[, 4]),   "\n", sep="")
  cat("                         Number of trees: ", length(x$weights),                     "\n",sep="")
  cat("                                minsplit: ", x$info$control$minsplit,               "\n", sep="")
  cat("                               minbucket: ", x$info$control$minbucket,              "\n", sep="")
  cat("                                 minprob: ", x$info$control$minprob,                "\n", sep="")
  cat("                                maxdepth: ", x$info$control$maxdepth,               "\n", sep="")
  cat("    No. of variables tried at each split: ", x$info$control$mtry,                   "\n", sep="")
  cat("                  Total no. of variables: ", length(attr(x$terms,'term.labels')),   "\n", sep="")
  cat("            Bootstrap type to grow trees: ", x$info$bootstrap,                      "\n",sep="")
  if (x$info$bootstrap != "by.user"){
    cat("           Resampling used to grow trees: ", samptype,                            "\n",sep="")
    if (samptype == "sampling without replacement"){
      cat("      Resampling rate used to grow trees: ", x$info$sampfrac,                   "\n",sep="")
    }
  }
  cat("                                Analysis: ", familyPretty,                 "\n", sep="")
  cat("                                  Family: ", familyOrg,                    "\n", sep="")
  
  cat("                          Splitting rule: ", "conditional inference framework",    "\n", sep="")
  
  
  #################################################################################
  ##
  ## end of grow mode
  ##
  #################################################################################
  
  
}

#' @export
print.ltrcrrf <- function(x) {
  # ## default printing
  # if (sum(inherits(x, c("ltrcrrf", "forest"), TRUE) == c(1, 2)) == 2) {
  #   print.default(x)
  #   return()
  # }
  
  ## check that the object is interpretable
  # if (sum(inherits(x, c("ltrcrrf", "grow"), TRUE) == c(1, 2)) != 2) {
  #   stop("This function only works for objects of class `(ltrcrrf, grow)'.")
  # }
  # grow.mode <- TRUE
  
  ## x will be processed if it's multivariate - therefore save some values from it
  familyPretty <- familypretty(x)
  familyOrg <- x$family
  yvar.dim <- ncol(x$yvar)
  
  
  ## ensure backward compatibility for nsplit
  if (is.null(x$nsplit)) {
    x$nsplit <- 0
  }
  
  if (x$forest$samptypeLTRC == "swor"){
    samptype <- "sampling without replacement"
  } else if (x$forest$samptypeLTRC == "swr"){
    samptype <- "sampling with replacement"
  }
  #################################################################################
  ##
  ## grow mode
  ##
  #################################################################################
  
  cat(" Number of (pseudo-subject) observations: ", x$n,                        "\n", sep="")
  
  cat("                      Number of subjects: ", length(unique(x$yvarLTRC$id)),           "\n",  sep="")
  cat("                        Number of deaths: ", sum(x$yvarLTRC[4]),   "\n", sep="")
  
  if (!is.null(x$imputed.indv)) {
    cat("                      Was data imputed: ", "yes",               "\n", sep="")
    #cat("                         Missingness: ",
    #    round(100*length(x$imputed.indv)/x$n,2), "%\n", sep="")
  }
  cat("                         Number of trees: ", x$ntree,                "\n",sep="")
  cat("               Forest terminal node size: ", x$nodesize,             "\n", sep="")
  cat("           Average no. of terminal nodes: ", mean(x$leaf.count),     "\n", sep="")
  cat("    No. of variables tried at each split: ", x$mtry,                 "\n", sep="")
  cat("                  Total no. of variables: ", length(x$xvar.names),   "\n", sep="")
  cat("            Bootstrap type to grow trees: ", x$forest$bootstrapLTRC,             "\n",sep="")
  if (x$forest$bootstrap != "by.user"){
    cat("           Resampling used to grow trees: ", samptype,            "          \n",sep="")
    if (samptype == "sampling without replacement"){
      cat("      Resampling rate used to grow trees: ", x$forest$sampfracLTRC,            "\n",sep="")
    }
  }
  cat("                                Analysis: ", "LTRCRRF",                  "\n", sep="")
  cat("                                  Family: ", familyOrg,                    "\n", sep="")
  if (x$nsplit > 0 & x$splitrule != "random") {
    cat("                          Splitting rule: ", paste(x$splitruleLTRC,"*random*"), "\n", sep="")
    cat("           Number of random split points: ", x$nsplit,  "\n", sep="")
  } else {
    cat("                          Splitting rule: ", x$splitruleLTRC,                   "\n", sep="")
  }
  
  #################################################################################
  ##
  ## end of grow mode
  ##
  #################################################################################
  
  
}
