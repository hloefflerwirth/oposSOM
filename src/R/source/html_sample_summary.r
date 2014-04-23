



outfile = file( paste( files.name, " - Results/Summary Sheets - Samples/0verview.html",sep=""), "w" )
            
  cat( "
    <html> <head> <TITLE>Sample Summary of ",files.name,"dataset</TITLE>
    </head> <body bgcolor=#FFFFFF >
    <style type=text/css> p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px }
    </style>

    <H1>Error Model</H1>       
    <TABLE BORDER=2, WIDTH=90%>
    <colgroup>
      <col width=\"50%\">
      <col width=\"50%\">
    </colgroup>", sep="", file = outfile)


  if( preferences$error.model == "replicates" )
  {
    cat( "
      <TR><TD>Error estimation model:</TD>
      <TD>Shrinkage: SD(replicates) vs. LPE</TD></TR>", sep="", file = outfile)
  } else
  if( preferences$error.model == "all.samples" )
  {
    cat( "
      <TR><TD>Error estimation model:</TD>
      <TD>LPE: SD(all genes) vs. <e>(all genes)</TD></TR>", sep="", file = outfile)
  } else
  if( preferences$error.model == "groups" )
  {
    cat( "
      <TR><TD>Error estimation model:</TD>
      <TD>Shrinkage: SD(categories) vs. LPE</TD></TR>", sep="", file = outfile)
  }


  cat( "
    <TR>
    <TD>LPE error plot for all samples:</TD>
    <TD><a href=\"../LPE/Sigma_LPE.pdf\" target=\"_blank\">PDF</a></TD>
    </TR></TABLE><br>




       
    <H1>Sample Summary Sheets</H1>
    
    <TABLE BORDER=2, WIDTH=90%>
    <TR><TD align=\"justify\">   
    For each sample a report sheet is created which summarizes the most relevant information using the global and local perspective.
    The global summary shows the ranked list of differentially expressed genes together with the associated significance characteristics for the whole sample, the ranked list of over- and underexpressed gene sets after GSZ-overexpression analysis and the respective p-value distributions.
    The local summary sheets present the analogous information for each single spot which is selected using the 98%-quantile criterion. The two maps in the left part of the sheet show the respective first level SOM and the selected spot, respectively.
    The full global and local lists can be downloaded in excel format.</TD>
    </TD></TR></TABLE><br>
       
    <TABLE BORDER=2, WIDTH=90%>
    <colgroup>
      <col width=\"22%\">
      <col width=\"18%\">
      <col width=\"12%\">
      <col width=\"12%\">
      <col width=\"36%\">
    </colgroup>
    <thead>
      <tr>
      <th>Sample name</th>
      <th>Category</th>
      <th>Summary Sheets</th>
      <th>Global Gene List</th>
      <th>Local Gene Lists</th>
      </tr>
    </thead>", sep="", file = outfile)


# <TD><a href=\"LPE/", colnames(indata)[m], ".bmp\" target=\"_blank\">BMP</a></TD>
  for( m in 1:ncol(indata) )
  {

    cat( "<TR>
      <TD>", colnames(indata)[m], "</TD>
      <TD>", group.labels[m], "</TD>
      <TD><a href=\"", colnames(indata)[m], ".pdf\" target=\"_blank\">PDF</a></TD>
      <TD><a href=\"../CSV Sheets/Gene Lists - Global/", colnames(indata)[m], ".csv\" >CSV</a></TD>

      <TD>", sep="", file = outfile)

  
    for( spot.i in 1:length( GS.infos.samples[[m]]$spots ) )
    {
      cat( "<a href=\"../CSV Sheets/Gene Lists - Local/", colnames(indata)[m], ".", spot.i, ".csv\" >CSV ",spot.i,"</a>&nbsp;&nbsp;&nbsp;", sep="", file = outfile)

    }

    cat( "</TD></TR>", sep="", file = outfile)

  }




  cat("  </table></body> </html> ", sep="", file = outfile)
close(outfile)






