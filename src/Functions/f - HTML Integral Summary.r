
	

outfile = file( paste( files.name, " - Results/Summary Sheets - Integral/0verview.html",sep=""), "w" )
          	
	cat( "
		<html> <head> <title>Spot Summary of ",files.name,"dataset</title>
		</head> <body bgcolor=#FFFFFF >
		<style type=text/css> p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px }
		</style>

	<h1>Integral Maps Summary Sheets</h1>
			 
			 
	<TABLE BORDER=2, WIDTH=90%>
	<TR><TD align=\"justify\">	 
	These analyses apply different criteria of spot selection such as overexpression, underexpression, maximum and minimum of metagene expression. Additionaly, mutual correlations between the metagenes and k-means clustering was applied to define spots of co-regulated metagenes.
	Gene set enrichment analysis provides the leading sets of each of the spots considered. Spot-related heatmaps characterize the expression profiles of the selected features in the series of samples.
	Single spot summary sheets provide detailed information on each of the spots such as the ranked list of samples which overexpress this feature in decreasing order according to the mean t-shrinkage statistics of the spot and the ranked list of the overrepresented gene sets together with the histogram of the respective p-value distribution.</TD>
	</TD></TR></TABLE><br>
		 
			 
	<table border=2, width=90%>
		<colgroup>
			<col width=\"25%\">
			<col width=\"20%\">
			<col width=\"55%\">
		</colgroup>
		<thead>
			<tr>
			<th>Type</th>
			<th>Summary Sheets</th>
			<th>Spot's CSV tables</th>
			</tr>
		</thead>


	<tr>
	<td>Overexpression</td>
	<td><a href=\"Overexpression.pdf\" target=\"_blank\">PDF</a></td>
	<td>&nbsp", sep="", file = outfile)

	for( m in 1:length(GS.infos.overexpression$spots) )
	{
		cat( "<a href=\"../CSV Sheets/Spot Lists/Sample-Overexpression ",names(GS.infos.overexpression$spots)[m],".csv\" target=\"_blank\">",names(GS.infos.overexpression$spots)[m],"</a>&nbsp;&nbsp;", sep="", file = outfile)
	}

	cat(
	"</td></tr>	


	<tr>
	<td>Underexpression</td>
	<td><a href=\"Underexpression.pdf\" target=\"_blank\">PDF</a></td>
	<td>&nbsp", sep="", file = outfile)

	for( m in 1:length(GS.infos.underexpression$spots) )
	{
		cat( "<a href=\"../CSV Sheets/Spot Lists/Sample-Underexpression ",names(GS.infos.underexpression$spots)[m],".csv\" target=\"_blank\">",names(GS.infos.underexpression$spots)[m],"</a>&nbsp;&nbsp;", sep="", file = outfile)
	}


# 	cat(
# 	"</td></tr>	
# 
# 
# 	<tr>
# 	<td>Metagene Maxima</td>
# 	<td><a href=\"Positive Metagenepeaks.pdf\" target=\"_blank\">PDF</a></td>
# 	<td>&nbsp", sep="", file = outfile)
# 
# 	for( m in 1:length(GS.infos.positivepeaks$spots) )
# 	{
# 		cat( "<a href=\"../CSV Sheets/Spot Lists/Metagene Maxima ",names(GS.infos.positivepeaks$spots)[m],".csv\" target=\"_blank\">",names(GS.infos.positivepeaks$spots)[m],"</a>&nbsp;&nbsp;", sep="", file = outfile)
# 	}
# 
# 	cat(
# 	"</td></tr>	
# 
# 
# 	<tr>
# 	<td>Metagene Minima</td>
# 	<td><a href=\"Negative Metagenepeaks.pdf\" target=\"_blank\">PDF</a></td>
# 	<td>&nbsp", sep="", file = outfile)
# 
# 	for( m in 1:length(GS.infos.negativepeaks$spots) )
# 	{
# 		cat( "<a href=\"../CSV Sheets/Spot Lists/Metagene Minima ",names(GS.infos.negativepeaks$spots)[m],".csv\" target=\"_blank\">",names(GS.infos.negativepeaks$spots)[m],"</a>&nbsp;&nbsp;", sep="", file = outfile)
# 	}

	cat(
	"</td></tr>	



	<tr>
	<td>Mutual Correlation</td>
	<td><a href=\"Correlation Cluster.pdf\" target=\"_blank\">PDF</a></td>
	<td>&nbsp", sep="", file = outfile)

	for( m in 1:length(GS.infos.correlation$spots) )
	{
		cat( "<a href=\"../CSV Sheets/Spot Lists/Correlation Cluster ",names(GS.infos.correlation$spots)[m],".csv\" target=\"_blank\">",names(GS.infos.correlation$spots)[m],"</a>&nbsp;&nbsp;", sep="", file = outfile)
	}

	cat(
	"</td></tr>

	
	
	<tr>
	<td>k-Means Clustering</td>
	<td><a href=\"K-Means Cluster.pdf\" target=\"_blank\">PDF</a></td>
	<td>&nbsp", sep="", file = outfile)

	for( m in 1:length(GS.infos.kmeans$spots) )
	{
		cat( "<a href=\"../CSV Sheets/Spot Lists/K-Means Cluster ",names(GS.infos.kmeans$spots)[m],".csv\" target=\"_blank\">",names(GS.infos.kmeans$spots)[m],"</a>&nbsp;&nbsp;", sep="", file = outfile)
	}


	cat(
	"</td></tr>

	
	
	<tr>
	<td>Group Overexpression</td>
	<td><a href=\"Group Overexpression.pdf\" target=\"_blank\">PDF</a></td>
	<td>&nbsp", sep="", file = outfile)

	for( m in 1:length(GS.infos.group.overexpression$spots) )
	{
		cat( "<a href=\"../CSV Sheets/Spot Lists/Group Overexpression ",names(GS.infos.group.overexpression$spots)[m],".csv\" target=\"_blank\">",names(GS.infos.group.overexpression$spots)[m],"</a>&nbsp;&nbsp;", sep="", file = outfile)
	}

	cat(
	"</td></tr></TABLE><br>



	
	<h1>Spot Report Sheets</h1>
			 
			 
	<TABLE BORDER=2, WIDTH=90%>
	<TR><TD align=\"justify\">	 
	Spot repots contain assignment of the spots to featuring samples and to respective sample groups. 
	Additionally, relative numbers of significant gene sets are given for each spot and gene set category.</TD>
	</TD></TR></TABLE><br>
		 
			 
	<table border=2, width=90%>
		<colgroup>
			<col width=\"50%\">
			<col width=\"50%\">
		</colgroup>
		<thead>
			<tr>
			<th>Type</th>
			<th>Summary Sheets</th>
			</tr>
		</thead>	
	
	<tr>
	<td>Overexpression</td>
	<td><a href=\"../3rd lvl Spot Analysis/Overexpression Spot Report.pdf\" target=\"_blank\">PDF</a></td>
	</tr>	

	<tr>
	<td>Underexpression</td>
	<td><a href=\"../3rd lvl Spot Analysis/Underexpression Spot Report.pdf\" target=\"_blank\">PDF</a></td>
	</tr>
	
	<tr>
	<td>K-Means</td>
	<td><a href=\"../3rd lvl Spot Analysis/K-Means Cluster Report.pdf\" target=\"_blank\">PDF</a></td>
	</tr>

	</TABLE>
	
	
	
	
	
	
	<h1>Spot Network Analysis</h1>
			 
			 
	<TABLE BORDER=2, WIDTH=90%>
	<TR><TD align=\"justify\">	 
	Networks of spot association are visualized as graphs. Both, correlation networks (connecting spots with r>0.5) and correlation spanning trees, are given for individual spots and also for patterns of spots occuring in the sample SOM images.</TD>
	</TD></TR></TABLE><br>
		 
			 
	<table border=2, width=90%>
		<colgroup>
			<col width=\"50%\">
			<col width=\"50%\">
		</colgroup>
		<thead>
			<tr>
			<th>Type</th>
			<th>Summary Sheets</th>
			</tr>
		</thead>	
	
	<tr>
	<td>Overexpression</td>
	<td><a href=\"../3rd lvl Spot Analysis/Overexpression Networks.pdf\" target=\"_blank\">PDF</a></td>
	</tr>	

	<tr>
	<td>Underexpression</td>
	<td><a href=\"../3rd lvl Spot Analysis/Underexpression Networks.pdf\" target=\"_blank\">PDF</a></td>
	</tr>

	<tr>
	<td>K-Means</td>
	<td><a href=\"../3rd lvl Spot Analysis/K-Means Cluster Networks.pdf\" target=\"_blank\">PDF</a></td>
	</tr>

	</TABLE>
	
	
	
	
	
	
	
	<h1>Chromosomal Enrichment</h1>
			 
			 
	<TABLE BORDER=2, WIDTH=90%>
	<TR><TD align=\"justify\">	 
	For each spot, enrichment of chromosomal postitions (chromosome/band) is visualized as overview heatmaps and individual chromosome plots.</TD>
	</TD></TR></TABLE><br>
		 
			 
	<table border=2, width=90%>
		<colgroup>
			<col width=\"50%\">
			<col width=\"50%\">
		</colgroup>
		<thead>
			<tr>
			<th>Type</th>
			<th>Summary Sheets</th>
			</tr>
		</thead>	
	
	<tr>
	<td>Overexpression</td>
	<td><a href=\"../3rd lvl Spot Analysis/Overexpression Chromosomal Enrichment.pdf\" target=\"_blank\">PDF</a></td>
	</tr>	

	<tr>
	<td>Underexpression</td>
	<td><a href=\"../3rd lvl Spot Analysis/Underexpression Chromosomal Enrichment.pdf\" target=\"_blank\">PDF</a></td>
	</tr>
	
	<tr>
	<td>K-Means</td>
	<td><a href=\"../3rd lvl Spot Analysis/K-Means Cluster Chromosomal Enrichment.pdf\" target=\"_blank\">PDF</a></td>
	</tr>

	</TABLE>	
	
	
	
	
	
	
	</body> </html> ", sep="", file = outfile)
close(outfile)











