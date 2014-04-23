# Installs CRAN packages if necessary
require.cran = function( package )
{
  if( length( find.package( package, quiet=T ) ) == 0 )
  {
    install.packages( package, repos="http://cran.r-project.org" )
  }
  suppressMessages( library( package, character.only=T, verbose=F ) )
}

# Installs Bioconductor packages if necessary
require.bioconductor = function( package )
{
  if( length( find.package( package, quiet=T ) ) == 0 )
  {
    source("http://bioconductor.org/biocLite.R")
    biocLite( package )
  }
  suppressMessages( library( package, character.only=T, verbose=F ) )
}
