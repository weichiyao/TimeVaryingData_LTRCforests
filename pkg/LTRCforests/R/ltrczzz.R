#' @importFrom parallel mclapply
#' @importFrom stats as.dist as.formula cutree dlnorm formula hclust lowess median model.matrix na.omit optim pgamma plnorm pnorm predict qnorm runif sd supsmu var wilcox.test
#' @importFrom utils installed.packages txtProgressBar setTxtProgressBar write.table tail
.onAttach <- function(libname, pkgname) {
  LTRCforests.version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                                  fields="Version")
  packageStartupMessage(paste("\n",
                              pkgname,
                              LTRCforests.version,
                              "\n",
                              "\n",
                              # "Type LTRCforests.news() to see new features, changes, and bug fixes.",
                              "\n",
                              "\n"))
}
