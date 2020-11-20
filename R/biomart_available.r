biomart.available <- function(env)
{
  mart <- try({ suppressMessages({ useMart(biomart=env$preferences$database.biomart, host=env$preferences$database.host, verbose=FALSE) }) }, silent=TRUE)

  return(is(mart,"Mart"))
}
