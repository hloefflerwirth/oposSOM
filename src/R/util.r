# outputs a message
util.cat <- function(..., suffix="\n") {
  message(paste(..., collapse=" "), suffix, appendLF=FALSE)
  flush.console()
}

# logs a message
util.log <- function(..., prefix, suffix="\n") {
  util.cat(c(paste("[", format(Sys.time()), "][", prefix, "]", sep=""), ...),
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
  x <- 1

  if (max <= len) {
    x <- round((cur * len) / max)
  } else {
    x <- which(floor(seq(0, max, length.out=len+1)) == floor(cur))[1] - 1
  }

  if (is.na(x)) {
    return()
  }

  util.info(paste("[",
                  paste(rep("#", x), collapse=""),
                  paste(rep(" ", len - x), collapse=""),
                  "]", sep=""), suffix="\r")
}

# terminate progress bar
util.progress.terminate <- function(len=48) {
  util.info(paste("[", paste(rep("#", len), collapse=""), "]", sep=""))
}

#
util.call <- function(fn, env) {
  environment(fn) <- env
  return(fn())
}

#
util.save <- function(env, filename) {
  return(tools:::makeLazyLoadDB(env, filename))
}

#
util.load <- function(filename) {
  env <- new.env()
  lazyLoad(filename, envir=env)
  return(env)
}
