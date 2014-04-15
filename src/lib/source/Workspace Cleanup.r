

remove.list = ls()[ which( !( ls() %in% c(

	## Functions
		
			
		"GeneSet.Fisher",
		"GeneSet.GSZ",
		
		"heatmap",
		"heatmap.wrap",
		
		"Get.Area",
		
		"Get.Running.Average",
		
		"get.beta.statistic",
		
		"require.cran",
		"require.bioconductor",
		
		"col.pix",
		"get.neighbors",
		
		"plot.set.list",
		"plot.set.list.chromosomes",
		"plot.set.list.reports",
		"plot.set.list.networks",
		"csv.set.list",
		

		
	## Primary objects
				
		"preferences",
		"files.name",
				
		"indata",
		"indata.original",
		"indata.mean.level",
		"indata.sample.mean",
		
		"group.labels",
		"group.colors",
		"unique.group.colors",
		
		"som.result",
		"metadata",
				
		"GS.infos.samples",
					
		"GS.infos.overexpression",
		"GS.infos.underexpression",
# 		"GS.infos.positivepeaks",
# 		"GS.infos.negativepeaks",		
		"GS.infos.correlation",
		"GS.infos.kmeans",
		"GS.infos.group.overexpression",
						
		
	## Suppoting data objects
				
		"loglog.metadata",
		"WAD.metadata",
		
		"group.metadata",
		"loglog.group.metadata",
		"WAD.group.metadata",
				
		"gene.names",
		"gene.ids",
		"gene.descriptions",
		"gene.positions",
		"gene.positions.list",
		"gene.positions.table",
		"genes.coordinates",
		"som.nodes",
		"unique.protein.ids",
		"gs.def.list",
		"gs.def.list.categories",	
				
		"sd.g.m",
		"t.g.m",
		"p.g.m",
		"fdr.g.m",
		"Fdr.g.m",
		"WAD.g.m",
		
		"n.0.m",
		"perc.DE.m",
					
		"t.m",
		"p.m",
		"fdr.m",
			
		"group.bootstrap.score",
		
		
	
	## Various data
		
		"LETTERS",
		"letters",
		
		"colramp",
														
		"supersom.custom",
		"supersom.20",
					
		"metagene.filter.list",
										
		"batch.t.g.m",
		"gs.null.list",
		"null.scores",
		"null.culdensity",
		
		"output.paths"
			
) ) ) ]			
			

rm( list = remove.list )

gc()

