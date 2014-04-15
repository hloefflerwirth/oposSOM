
 

if( preferences$geneset.analysis )
{



	outfile = file( paste( files.name, " - Results/Geneset Analysis/0verview.html",sep=""), "w" )
          	
	cat( "
		<html> <head> <TITLE>Geneset Analysis Summary of",files.name,"dataset</TITLE>
		</head> <body bgcolor=#FFFFFF >
		<style type=text/css> p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px }
		</style>

			 
		<H1>General Information</H1>
		<TABLE BORDER=2, WIDTH=90%>
		<colgroup>
			<col width=\"50%\">
			<col width=\"50%\">
		</colgroup>

		<TR>
		<TD>Number of Genesets:</TD>
		<TD>", length(gs.def.list), "</TD>
		</TR>

		<TR>
		<TD>Categories:</TD>
		<TD>", paste( paste( names( table( gs.def.list.categories ) ), table( gs.def.list.categories ), sep=" (" ), collapse=") , " ) , ")</TD>
		</TR>

		<TR>
		<TD>Table of all GSZ scores:</TD>
		<TD><a href=\"../CSV sheets/Sample GSZ scores.csv\" target=\"_blank\">CSV</a></TD>
		</TR>
			 
	 	<TR>
		<TD>GSZ/Fisher analysis: Heatmaps & p-value histograms:</TD>
	 	<TD><a href=\"0verview Heatmaps.pdf\" target=\"_blank\">PDF</a></TD>
	 	</TR>			 

		</TABLE><br>", sep="", file = outfile)


	
	cat( "
			 <H1>Quick Links</H1>
			 Go to category: ", sep="", file = outfile)
	
	for( i in names( table( gs.def.list.categories ) ) )
	{
		cat( "<a href=\"#", i, "\">", i, "</a>&nbsp;&nbsp;", sep="", file = outfile )
	}

	

	cat( "<br><br><H1>Gene Sets</H1>
	
		<TABLE BORDER=2, WIDTH=90%>
		<TR><TD align=\"justify\">	 
		Enrichment profiles of individual predefined gene sets are shown as bar plots across all samples. Additionally the log FC-expression profiles of the leading metagenes are shown.
		Further, members of each gene set are shown as population maps and listed in Excel-files.</TD>
		</TD></TR></TABLE><br><br>", sep="", file = outfile)	

	
	

	for( i in names( table( gs.def.list.categories ) ) )
	{
		category.gs.list = gs.def.list[ which( gs.def.list.categories == i ) ]


		cat( "<b><u><a name=\"", i, "\">Category ", i, "</a></u></b><br><br>", sep="", file = outfile)

		cat( "<TABLE BORDER=2, WIDTH=90%>
			<colgroup><col width=\"55%\"><col width=\"9%\"><col width=\"12%\"><col width=\"12%\"><col width=\"12%\"></colgroup>
			<thead><tr>
				<th>Geneset name</th>
				<th>Category</th>
				<th>Profile</th>
				<th>Population Map</th>
				<th>Members</th>
			</tr>	</thead>", sep="", file = outfile)

		for( ii in 1:length(category.gs.list)  )
		{
			cat( "<TR>
				<TD>", names(category.gs.list)[ii], "</TD>
				<TD>", sapply(category.gs.list,function(x){x$Type})[ii], "</TD>
				<TD><a href=\"", substring( make.names( names(category.gs.list)[ii]), 1, 100 ), " profile.pdf\" target=\"_blank\">PDF</a></TD>
				<TD><a href=\"", substring( make.names( names(category.gs.list)[ii]), 1, 100 ), " map.pdf\" target=\"_blank\">PDF</a></TD>
				<TD><a href=\"", substring( make.names( names(category.gs.list)[ii]), 1, 100 ), ".csv\" >CSV</a></TD>
				</TR>", sep="", file = outfile)
		}

		cat( "</TABLE><br><br>", sep="", file = outfile)

	}





	cat("	</body> </html> ", sep="", file = outfile)
	close(outfile)
	
}


 
