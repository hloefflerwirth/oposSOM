biomart.available <- function()
{
  mart <- try({ suppressMessages({ useMart(biomart=preferences$database.biomart, host=preferences$database.host, verbose=FALSE) }) }, silent=TRUE)

  return(class(mart)=="Mart")
}
