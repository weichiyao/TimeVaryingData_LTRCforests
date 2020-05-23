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
