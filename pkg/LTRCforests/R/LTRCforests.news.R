LTRCforests.news <- function(...) {
  newsfile <- file.path(system.file(package="LTRCforests"), "NEWS")
  file.show(newsfile)
}

