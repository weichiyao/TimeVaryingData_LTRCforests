#' Print Summary Output of a ltrccf object
#'
#' Print summary output after a LTRCCF model is built.
#' This is the default print method for objects in the class of \code{\link{ltrccf}}.
#'
#' @param x an object of class \code{(ltrccf, grow)}.
#'
#' @importFrom survival Surv
#' @examples
#'
#' library(survival)
#' Formula = Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' # Built a LTRCCF forest on the time-varying data by specifying id, with mtry specified:
#' LTRCCFobj = ltrccf(formula = Formula, data = pbcsample, id = ID, mtry = 3, ntree = 50L)
#' print(LTRCCFobj)
#'
#' # Built a LTRCCF forest on the time-invariant data, with resampling, with mtry specified:
#' LTRCCFobj = ltrccf(formula = Formula, data = pbcsample, samptype = "swr",
#'                    mtry = 3, ntree = 50L)
#' print(LTRCCFobj)
#' @seealso \code{\link{ltrccf}}
#' @export



print.ltrccf <- function(x) {
  ## check that the object is interpretable
  if (any(inherits(x, c("ltrccf", "grow"), TRUE) ) == 0 ) {
    stop("This function only works for objects (of subclass) of class `(ltrccf, grow).")
  }
  grow.mode <- TRUE

  ## x will be processed if it's multivariate - therefore save some values from it
  familyPretty <- "LTRCCF"
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
