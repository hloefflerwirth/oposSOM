util.cat <- function(...) {
  cat(c(...), "\n")
  flush.console()
}

util.log <- function(..., prefix) {
  util.cat(paste("[", format(Sys.time()), "][", prefix, "]", sep=""), c(...))
}

util.info <- function(...) {
  util.log(c(...), prefix="INFO")
}

util.warn <- function(...) {
  util.log(c(...), prefix="WARN")
}

util.fatal <- function(...) {
  util.log(c(...), prefix="FATAL")
}
