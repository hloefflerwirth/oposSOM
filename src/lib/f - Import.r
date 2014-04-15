# Import libraries
source( "lib/quantile_normalization.r" )
source( "lib/geneset_analysis.r" )
source( "lib/get_running_average.r" )
source( "lib/get_area.r" )
source( "lib/get_statistic.r" )
source( "lib/heatmaps.r" )

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

# Load external packages
require.cran( "som" )
require.cran( "fastICA" )
require.cran( "scatterplot3d" )
require.cran( "pixmap" )
require.cran( "fdrtool" )
require.cran( "igraph" )
require.cran( "ape" )
require.cran( "KernSmooth" )
require.cran( "parallel" )
require.cran( "foreach" )

if( Sys.info()["sysname"] == "Windows" ) require.cran("doSNOW") else require.cran("doMC")
