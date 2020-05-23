print.ltrcrsf <- function(x, outcome.target = NULL, ...) {
  ## default printing
  if (sum(inherits(x, c("ltrcrsf", "forest"), TRUE) == c(1, 2)) == 2) {
    print.default(x)
    return()
  }

  ## check that the object is interpretable
  if (sum(inherits(x, c("ltrcrsf", "grow"), TRUE) == c(1, 2)) != 2 ) {
    stop("This function only works for objects of class `(ltrcrsf, grow).")
  }
  grow.mode <- TRUE

  ## x will be processed if it's multivariate - therefore save some values from it
  familyPretty <- family.pretty(x)
  familyOrg <- x$family
  yvar.dim <- ncol(x$yvar)

  ## survival: event frequencies
  if (grepl("surv", x$family)) {
    event <- get.event.info(x)$event
    n.event <- 1
    if (!is.null(event)) {
      n.event <- length(unique(event))
      if (length(event) > 0) {
        event.freq <- paste(tapply(event, event, length), collapse = ", ")
      }
        else {
          event.freq <- 0
        }
    }
  }

  ## error rates
  if (!is.null(x$err.rate)) {
    err.rate <- cbind(x$err.rate)
    if (grepl("surv", x$family)) {
      err.rate <- paste(round(100 * err.rate[nrow(err.rate), ], 2), "%", collapse=", ", sep = "")
    }
    else if (x$family == "class") {
      brierS <- round(100 * brierS, 2)
      aucS <- round(100 * aucS, 2)
      overall.err.rate <- paste(round(100 * err.rate[nrow(err.rate), 1], 2), "%", sep = "")
      ## rfq related adjustments
      if (!is.null(gmeanS)) {
        err.rate <- round(err.rate[nrow(err.rate), 1], 2)
      }
      else {
        err.rate <- paste(round(err.rate[nrow(err.rate), ], 2), collapse=", ", sep = "")
      }
    }
    else if (x$family == "regr") {
      per.var <- round(100 * (1 - err.rate[nrow(err.rate), ] / var(x$yvar, na.rm = TRUE)), 2)
      err.rate <- round(err.rate[nrow(err.rate), ], 2)
    }
    else {
      err.rate <- NULL
    }
  }
  else {
    err.rate <- NULL
  }
  ## ensure backward compatibility for nsplit
  if (is.null(x$nsplit)) {
    x$nsplit <- 0
  }
  #################################################################################
  ##
  ## grow mode
  ##
  #################################################################################

  cat("                         Sample size: ", x$n,                 "\n", sep="")
  if (grepl("surv", x$family)) {
    if (n.event > 1) {
      cat("                    Number of events: ", event.freq,  "\n", sep="")
    }
      else {
        cat("                    Number of deaths: ", x$ndead,   "\n", sep="")
      }
  }

  if (!is.null(x$imputed.indv)) {
    cat("                    Was data imputed: ", "yes",               "\n", sep="")
    #cat("                         Missingness: ",
    #    round(100*length(x$imputed.indv)/x$n,2), "%\n", sep="")
  }
  cat("                     Number of trees: ", x$ntree,                "\n",sep="")
  cat("           Forest terminal node size: ", x$nodesize,             "\n", sep="")
  cat("       Average no. of terminal nodes: ", mean(x$leaf.count),     "\n", sep="")
  cat("No. of variables tried at each split: ", x$mtry,                 "\n", sep="")
  cat("              Total no. of variables: ", length(x$xvar.names),   "\n", sep="")
  if (!x$univariate) {
    cat("              Total no. of responses: ", yvar.dim,   "\n", sep="")
    cat("         User has requested response: ", outcome.target,         "\n", sep="")
  }
  cat("        Bootstrap type to grow trees: ", x$bootstrap,                  "\n",sep="")
  cat("       Resampling used to grow trees: ", x$forest$samptype,            "\n",sep="")
  cat("    Resample size used to grow trees: ", round(x$forest$sampsize(x$n)),"\n",sep="")
  cat("                            Analysis: ", familyPretty,                 "\n", sep="")
  cat("                              Family: ", familyOrg,                    "\n", sep="")
  if (x$nsplit > 0 & x$splitrule != "random") {
    cat("                      Splitting rule: ", paste(x$splitrule,"*random*"),"\n", sep="")
    cat("       Number of random split points: ", x$nsplit                   ,  "\n", sep="")
  }
    else {
      cat("                      Splitting rule: ", x$splitrule,          "\n", sep="")
    }
  if (!is.null(err.rate)) {
    if (x$family == "regr") {
      cat("                % variance explained: ", per.var, "\n", sep="")
    }
    if (x$family == "class" && !is.null(brierS)) {
      cat("              Normalized brier score:", brierS, "\n")
    }
    if (x$family == "class" && !is.null(aucS)) {
      cat("                                 AUC:", aucS, "\n")
    }
    if (x$family == "class" && !is.null(gmeanS)) {
      cat("                              G-mean: ", gmeanS,            "\n", sep="")
      cat("                    Imbalanced ratio: ", iratio, "\n", sep="")
    }
    cat("                          Error rate: ", err.rate,            "\n\n", sep="")
  }
  #################################################################################
  ##
  ## end of grow mode
  ##
  #################################################################################




  if (sf.flag) {
    message(sf.message)
  }
}
