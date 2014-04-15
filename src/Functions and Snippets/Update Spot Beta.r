

source("lib/get_statistic.r")
require.cran( "fdrtool" )


i=1
for( i in 1:length(GS.infos.overexpression$spots) )		
	GS.infos.overexpression$spots[[i]]$beta.statistic = 
#		get.beta.statistic( set.data=metadata[ GS.infos.overexpression$spots[[ i ]]$metagenes,,drop=F ] )$beta.score				
		get.beta.statistic( set.data=metadata[ GS.infos.overexpression$spots[[ i ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.overexpression$spots[[ i ]]$metagenes, ]$nobs )				

for( i in 1:length(GS.infos.underexpression$spots) )		
	GS.infos.underexpression$spots[[i]]$beta.statistic = 
	get.beta.statistic( set.data=metadata[ GS.infos.underexpression$spots[[ i ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.underexpression$spots[[ i ]]$metagenes, ]$nobs )			


for( i in 1:length(GS.infos.positivepeaks$spots) )		
	GS.infos.positivepeaks$spots[[i]]$beta.statistic = 
	get.beta.statistic( set.data=metadata[ GS.infos.positivepeaks$spots[[ i ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.positivepeaks$spots[[ i ]]$metagenes, ]$nobs )	

for( i in 1:length(GS.infos.negativepeaks$spots) )		
	GS.infos.negativepeaks$spots[[i]]$beta.statistic = 
	get.beta.statistic( set.data=metadata[ GS.infos.negativepeaks$spots[[ i ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.negativepeaks$spots[[ i ]]$metagenes, ]$nobs )


for( i in 1:length(GS.infos.correlation$spots) )		
	GS.infos.correlation$spots[[i]]$beta.statistic = 
	get.beta.statistic( set.data=metadata[ GS.infos.correlation$spots[[ i ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.correlation$spots[[ i ]]$metagenes, ]$nobs )

for( i in 1:length(GS.infos.kmeans$spots) )		
	GS.infos.kmeans$spots[[i]]$beta.statistic = 
	get.beta.statistic( set.data=metadata[ GS.infos.kmeans$spots[[ i ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.kmeans$spots[[ i ]]$metagenes, ]$nobs )	





dir.create( paste( files.name, "- Results" ), showWarnings=F )
dir.create( paste( files.name, "- Results/CSV Sheets" ), showWarnings=F )

source("lib/source/Summary Sheets Integral.r")










