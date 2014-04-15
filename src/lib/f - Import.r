

source( "lib/l - Quantile Normalization.r" )
source( "lib/l - GeneSet Analysis.r" )
source( "lib/l - Get.Running.Average.r" )
source( "lib/l - Get.Area.r" )                                                 
source( "lib/l - Get Statistic.r" )                                                 	
source( "lib/l - Heatmaps.r" )    





require.cran = function( package )
{
	if( length( find.package( package, quiet=T ) ) == 0 )
	{
		install.packages( package, repos="http://cran.r-project.org" )
	}		
	suppressMessages( library( package, character.only=T, verbose=F ) )
}

require.bioconductor = function( package )
{
	if( length( find.package( package, quiet=T ) ) == 0 )
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite( package )
	}
	suppressMessages( library( package, character.only=T, verbose=F ) ) 
}


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


