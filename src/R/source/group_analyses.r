

perform.group.analyses = function()
{
	cat( "Group-centered Analyses\n\n" ); flush.console()
	
	dir.create( paste( files.name, "- Results/Summary Sheets - Groups" ), showWarnings=F )		
	dir.create( paste( files.name, "- Results/Summary Sheets - Groups/CSV Sheets" ), showWarnings=F )

	if( preferences$geneset.analysis ) 
	{
		dir.create( paste( files.name, "- Results/Summary Sheets - Groups/Geneset Analysis" ), showWarnings=F )
		source("R/source/group_specific_genesets.r")
	}	
	
	
	
	source("R/source/summary_sheets_groups.r")	

	
	
	metadata = do.call( cbind, by( t(metadata), group.labels, colMeans )[unique(group.labels)] )
			
	
	colnames(indata.original) = group.labels[ colnames(indata.original) ]
	
	indata = do.call( cbind, by( t(indata.original), colnames(indata.original), colMeans )[unique(group.labels)] )
		
	indata.mean.level = rowMeans( indata )
	
	if( preferences$feature.mean.normalization )
	{		
		indata = indata - indata.mean.level	
	}

	
	group.colors = group.colors[ match( colnames(indata), group.labels) ]	
	group.labels = group.labels[ match( colnames(indata), group.labels) ]
	names(group.labels) = group.labels
	names(group.colors) = group.labels
	
	output.paths = c( "LPE" = "",
										"CSV" = paste( files.name, "- Results/Summary Sheets - Groups/CSV Sheets" ),
										"Summary Sheets Samples"= paste( files.name, "- Results/Summary Sheets - Groups/Reports" ) )



	preferences = preferences
	preferences$error.model = "replicates"
	
	capture.output({ 
		source("R/source/calc_statistics.r", local=T);		

		source("R/source/detect_spots_samples.r", local=T);

		if( preferences$geneset.analysis ) 
		{
			source("R/source/geneset_statistic_samples.r", local=T);
		}

		
		source("R/source/gene_lists.r", local=T);

		source("R/source/summary_sheets_samples.r", local=T);
	})
	
}




if( length(unique( group.labels ) ) > 1 )	
{
	suppressWarnings( perform.group.analyses() )


	if( length(unique( group.labels ) ) < 8 )	
	{
		differences.list = apply( combn( unique(group.labels), 2 ), 2, function(x)
		{
			list( which(group.labels==x[1]), which(group.labels==x[2]) ) 
		} )
		names(differences.list)  = apply( combn( unique(group.labels), 2 ), 2, paste, collapse=" vs " )
	
		preferences$differences.list = c( preferences$differences.list, differences.list )
	}

}	
	
	



