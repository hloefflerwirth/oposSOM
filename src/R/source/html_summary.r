


outfile = file( paste( files.name, " - Results/Summary.html",sep=""), "w" )
            
cat( "
  <html> <head> <TITLE>Summary of ",files.name," dataset</TITLE>
  </head> <body bgcolor=#FFFFFF >
  <style type=text/css> p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px }
  </style> ", file = outfile)




cat( "<H1>General Information</H1>
  <TABLE BORDER=2, WIDTH=90%>
  <colgroup>
    <col width=\"50%\">
    <col width=\"50%\">
  </colgroup>

  <TR>
  <TD>Number of samples:</TD>
  <TD>", ncol(indata), " (all replicates: ", ncol(indata.original), ")</TD>
  </TR>

  <TR>
  <TD>Number of categories:</TD>
  <TD>", length(unique(group.labels)), "</TD>
  </TR>

  <TR>
  <TD>Number of genes:</TD>
  <TD>", nrow(indata), "</TD>
  </TR>

  <TR>
  <TD>Analysis finished:</TD>
  <TD>", format(Sys.time(), "%a %b %d %X %Y %Z"), "</TD>
  </TR>

  </TABLE><br>", sep="", file = outfile)




cat( "<TABLE BORDER=2, WIDTH=90%>
  <colgroup>
    <col width=\"50%\">
    <col width=\"50%\">
  </colgroup>

  <TR>
  <TD>Dimension 1st level SOM:</TD>
  <TD>", preferences$dim.som1, "x", preferences$dim.som1, "</TD>
  </TR>

  <TR>
  <TD>Dimension 2nd level SOM:</TD>
  <TD>", preferences$dim.som2, "x", preferences$dim.som2, "</TD>
  </TR>

  </TABLE>", sep="", file = outfile)







cat( "<H1>Results</H1>

  <TABLE BORDER=2 , WIDTH=90%>
  <colgroup>
    <col width=\"50%\">
    <col width=\"50%\">
  </colgroup>
  <thead><tr>
    <th colspan=2, BGCOLOR=\"#99CCFF\" >1st Level SOM Analysis</th>
  </tr>  </thead>

  
  <TR>
  <TD rowspan=5 >These reports show the collection of first level SOM images of all tissue samples. Alternative SOM images are shown with different color scales and also as rank SOM using FC-, WAD-, and Shrinkage-t-scores.
     Further, supporting maps, entropy and topology profiles provide supplementary information about the 1st level SOM.</TD>

  <TD><a href=\"Expression Portraits.pdf\" target=\"_blank\">
    1st level SOM expression images</a>
  &nbsp;&nbsp;&nbsp;
  <a href=\"Expression Portraits alternative.pdf\" target=\"_blank\">
    alternative color scales</a>
  &nbsp;&nbsp;&nbsp;
  <a href=\"Rank Maps.pdf\" target=\"_blank\">
    rank profiles</a></TD>
  </TR>


  <TR>
   <TD><a href=\"Summary Sheets - Groups/Expression Portraits Groups.pdf\" target=\"_blank\">
     1st level SOM expression images - group specific</a></TD>
   </TR>
  <TR>

  <TR>
   <TD><a href=\"Supporting Maps&Profiles/Supporting Maps.pdf\" target=\"_blank\">
     Supporting Maps</a></TD>
   </TR>
  <TR>
     
   <TD><a href=\"Supporting Maps&Profiles/Entropy Profiles.pdf\" target=\"_blank\">
     Entropy Profiles</a>
   &nbsp;&nbsp;&nbsp;
   <a href=\"Supporting Maps&Profiles/Topology Profiles.pdf\" target=\"_blank\">
     Topology Profiles</a></TD>
   </TR>
 



     
  </TABLE><br>







  <TABLE BORDER=2, WIDTH=90%>
  <colgroup>
    <col width=\"50%\">
    <col width=\"50%\">
  </colgroup>
  <thead><tr>
    <th colspan=5, BGCOLOR=\"#99CCFF\" >Sample Summaries</th>
  </tr>  </thead>

  
  <TR>
  <TD rowspan=5 >For each sample a report sheet is created which summarizes the most relevant information.</TD>

  <TD><a href=\"Summary Sheets - Samples/0verview.html\" target=\"_blank\">
    browse sample report sheets</a></TD>
  </TR>

  </TABLE><br>", sep="", file = outfile)

     

     



  if( preferences$geneset.analysis )
  {
     
    cat("
    <TABLE BORDER=2, WIDTH=90%>
    <colgroup>
      <col width=\"50%\">
      <col width=\"50%\">
    </colgroup>
    <thead><tr>
      <th colspan=5, BGCOLOR=\"#99CCFF\" >Gene Set Enrichment Analysis</th>
    </tr>  </thead>
  
    
    <TR>
    <TD rowspan=5 >Enrichment of predefined gene sets in the samples is calculated and visualized as heatmaps, profile plots and population maps.</TD>
  
    <TD><a href=\"Geneset Analysis/0verview.html\" target=\"_blank\">
      browse geneset analysis results</a></TD>
  
    </TABLE><br>", sep="", file = outfile)
        
     
  }     
     
     
     
  cat("
  <br><br>
  <TABLE BORDER=2, WIDTH=90%>
  <colgroup>
    <col width=\"50%\">
    <col width=\"50%\">
  </colgroup>
  <thead><tr>
    <th colspan=4, BGCOLOR=\"#99CCFF\" >2nd Level Metagene Analysis</th>
  </tr>  </thead>

  
  <TR>
  <TD rowspan=5 >Several agglomerative methods based either on distance or on correlation metrics are applied to the samples using filtered subsets of metagenes.
      The reports show cluster heatmaps and dendrograms, pairwise correlation maps, correlation spanning trees and ICA results.<br>
      </TD>


   <TD><a href=\"2nd lvl Metagene Analysis/2nd lvl SOM.pdf\" target=\"_blank\">
    Second level SOM</a></TD>
   </TR>
 
  <TR>  
  <TD><a href=\"2nd lvl Metagene Analysis/Similarity Analysis.pdf\" target=\"_blank\">
    Similarity based methods ( Neighbor Joining & Hierarchical Clustering )</a></TD>
  </TR>

  <TR>
  <TD><a href=\"2nd lvl Metagene Analysis/Correlation Analysis.pdf\" target=\"_blank\">
    Correlation based methods ( MST, Correlation Graphs, PCM )</a></TD>
  </TR>

  <TR>
  <TD><a href=\"2nd lvl Metagene Analysis/Component Analysis.pdf\" target=\"_blank\">
    Component based methods ( 2d-ICA, 3d-ICA )</a></TD>
  </TR>


  </TABLE><br>


     
  <br><br>
  <TABLE BORDER=2, WIDTH=90%>
  <colgroup>
    <col width=\"50%\">
    <col width=\"50%\">
  </colgroup>
  <thead><tr>
    <th colspan=5, BGCOLOR=\"#99CCFF\" >3rd Level Spot Analysis</th>
  </tr>  </thead>

  
  <TR>
  <TD rowspan=6 >Different criteria of spot selection such as overexpression or mutual correlations between the metagenes where applied. 
     Summary sheets of integral SOM maps and subsequent analyses base on this third level information.</TD>

  <TD><a href=\"Summary Sheets - Integral/0verview.html\" target=\"_blank\">
    browse spot report sheets</a></TD>
  </TR>


  </TABLE><br>", sep="", file = outfile)






cat("  </body> </html> ", sep="", file = outfile)
close(outfile)




