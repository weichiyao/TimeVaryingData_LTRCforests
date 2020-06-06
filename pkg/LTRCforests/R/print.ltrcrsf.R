#' Print Summary Output of a ltrcrsf object
#'
#' Print summary output after a LTRCRSF model is built.
#' This is the default print method for objects in the class of \code{\link{ltrcrsf}}.
#'
#' @param x an object of class \code{(ltrcrsf, grow)}.
#'
#' @export
#' @examples
#' library(survival)
#' Formula = Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
#' # Built a LTRCRSF forest (based on bootstrapping subjects without replacement)
#' # on the time-varying data by specifying id, with mtry specified:
#' LTRCRSFobj = ltrcrsf(formula = Formula, data = pbcsample, id = ID, mtry = 3, ntree = 50L)
#' print(LTRCRSFobj)
#'
#' # Built a LTRCRSF forest (sampling with replacement) on the time-invariant data,
#' # with mtry specified:
#' LTRCRSFobj = ltrcrsf(formula = Formula, data = pbcsample, samptype = "swr", mtry = 3,
#'                      ntree = 50L)
#' print(LTRCRSFobj)
#' @seealso \code{\link{ltrcrsf}}

print.ltrcrsf <- function(x) {
  # ## default printing
  # if (sum(inherits(x, c("ltrcrsf", "forest"), TRUE) == c(1, 2)) == 2) {
  #   print.default(x)
  #   return()
  # }

  ## check that the object is interpretable
  if (sum(inherits(x, c("ltrcrsf", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(ltrcrsf, grow)'.")
  }
  grow.mode <- TRUE

  ## x will be processed if it's multivariate - therefore save some values from it
  familyPretty <- family.pretty(x)
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
  cat("                                Analysis: ", "LTRCRSF",                  "\n", sep="")
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
