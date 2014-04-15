

perform.group.analyses = function()
{
	cat( "Group-centered Analyses\n\n" ); flush.console()
	
	dir.create( paste( files.name, "- Results/Summary Sheets - Groups" ), showWarnings=F )		
	dir.create( paste( files.name, "- Results/Summary Sheets - Groups/CSV Sheets" ), showWarnings=F )

	if( preferences$geneset.analysis ) 
	{
		dir.create( paste( files.name, "- Results/Summary Sheets - Groups/Geneset Analysis" ), showWarnings=F )
		source("lib/f - Group Specific Genesets.r")
	}	
	
	
	
	source("lib/f - Summary Sheets Groups.r")	

	
	
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
		source("lib/f - Calc Statistics.r", local=T);		

		source("lib/f - Detect Spots Samples.r", local=T);

		if( preferences$geneset.analysis ) 
		{
			source("lib/f - Geneset Statistic Samples.r", local=T);
		}

		
		source("lib/f - Gene Lists.r", local=T);

		source("lib/f - Summary Sheets Samples.r", local=T);
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
	
	



