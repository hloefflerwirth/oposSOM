




	### init parallel computing ###

	try({ stopCluster(cl) }, silent=T )

	cl <- makeCluster( preferences$max.parallel.cores ) 
	

	


	### perform GS analysis ###

	
	for( i in 1:length( GS.infos.overexpression$spots ) )
	{
		geneset.ids = unique( na.omit( gene.ids[ GS.infos.overexpression$spots[[i]]$genes ] ) )
		GS.infos.overexpression$spots[[i]]$Fisher.p = GeneSet.Fisher( geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl )
	}
	cat( "#" );	flush.console()

	for( i in 1:length( GS.infos.underexpression$spots ) )
	{
		geneset.ids = unique( na.omit( gene.ids[ GS.infos.underexpression$spots[[i]]$genes ] ) )
		GS.infos.underexpression$spots[[i]]$Fisher.p = GeneSet.Fisher( geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl )
	}
	cat( "#" );	flush.console()
	
# 	for( i in 1:length( GS.infos.positivepeaks$spots ) )
# 	{
# 		geneset.ids = unique( na.omit( gene.ids[ GS.infos.positivepeaks$spots[[i]]$genes ] ) )
# 		GS.infos.positivepeaks$spots[[i]]$Fisher.p = GeneSet.Fisher( geneset.ids, unique.protein.ids, gs.def.list, sort=T )
# 	}
# 	
# 	for( i in 1:length( GS.infos.negativepeaks$spots ) )
# 	{
# 		geneset.ids = unique( na.omit( gene.ids[ GS.infos.negativepeaks$spots[[i]]$genes ] ) )
# 		GS.infos.negativepeaks$spots[[i]]$Fisher.p = GeneSet.Fisher( geneset.ids, unique.protein.ids, gs.def.list, sort=T )
# 	}
# 	cat( "#" );	flush.console()	

	if( length(GS.infos.correlation$spots) > 0 )
	for( i in 1:length( GS.infos.correlation$spots ) )
	{
		geneset.ids = unique( na.omit( gene.ids[ GS.infos.correlation$spots[[i]]$genes ] ) )
		GS.infos.correlation$spots[[i]]$Fisher.p = GeneSet.Fisher( geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl )
	}
	cat( "#" );	flush.console()

	for( i in 1:length( GS.infos.kmeans$spots ) )
	{
		geneset.ids = unique( na.omit( gene.ids[ GS.infos.kmeans$spots[[i]]$genes ] ) )
		GS.infos.kmeans$spots[[i]]$Fisher.p = GeneSet.Fisher( geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl )
	}
	cat( "#" );	flush.console()

	for( i in 1:length( GS.infos.group.overexpression$spots ) )
	{
		geneset.ids = unique( na.omit( gene.ids[ GS.infos.group.overexpression$spots[[i]]$genes ] ) )
		GS.infos.group.overexpression$spots[[i]]$Fisher.p = GeneSet.Fisher( geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl )
	}
	cat( "#" );	flush.console()
	
	
	cat("|\n");	flush.console()







	
	
	
	### stop parallel computing ###
	
	try({ stopCluster(cl) }, silent=T )
	
