# outputs a message
util.cat <- function(..., suffix="\n") {
  cat(c(...), suffix)
  flush.console()
}

# logs a message
util.log <- function(..., prefix, suffix="\n") {
  util.cat(paste("[", format(Sys.time()), "][", prefix, "]", sep=""),
           c(...),
           suffix=suffix)
}

# logs an info
util.info <- function(..., suffix="\n") {
  util.log(c(...), prefix="INFO", suffix=suffix)
}

# logs a warning
util.warn <- function(..., suffix="\n") {
  util.log(c(...), prefix="WARN", suffix=suffix)
}

# logs an fatal error
util.fatal <- function(..., suffix="\n") {
  util.log(c(...), prefix="FATAL", suffix=suffix)
}

# outputs a progress bar
util.progress <- function(cur, max, len=48) {
  # TODO: only cat() on change
  x <- round((cur * len) / max, 0)

  util.info(paste("[",
                  paste(rep("#", x), collapse=""),
                  paste(rep(" ", len - x), collapse=""),
                  "]", sep=""), suffix="\r")
}
