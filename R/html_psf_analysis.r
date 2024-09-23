pipeline.htmlPsfAnalysis <- function(env)
{
  dirname <- "PSF Analysis"

  if (!file.exists(dirname))
  {
    return()
  }
  
  data("kegg.collection")
  folders <- dir( "PSF Analysis" )
  
  if(length(folders)>0)
  for( x in folders )
  {
    
    filename <- file.path("PSF Analysis", x, "0verview.html" )
    util.info("Writing:", filename)
    outfile <- file(filename, "w")
  
    cat("<!DOCTYPE html>
    <html>
      <head>
        <title>PSF Analysis Summary of ", env$files.name, " dataset</title>
        <style>
          body {
            margin: 0;
            padding: 0;
            color: #333;
            background-color: #fff;
            font: normal normal normal 14px/180% sans-serif;
          }
    
          #wrapper {
            width: 90%;
            min-width: 400px;
            max-width: 800px;
            margin: 20px auto;
          }
    
          h1, h2 {
            clear: both;
            margin: 30px 0 0 0;
            line-height: 210%;
            border-bottom: 1px solid #eee;
          }
    
          dl {
            line-height: 180%;
          }
    
          dl dt {
            width: 50%;
            float: left;
            color: #111;
          }
    
          dl dt:after {
            content: ':';
          }
    
          dl dd {
            margin-left: 50%;
          }
    
          a {
            color: #4183C4;
            text-decoration: none;
          }
    
          table {
            width: 100%;
            margin: 24px 0;
            border-collapse: collapse;
          }
    
          table th,
          table td {
            text-align: left;
            padding: 4px 6px;
          }
    
          table thead tr,
          table tbody tr:nth-child(2n) {
            background-color: #f0f0f0;
          }
    
          table td a {
            display: block;
          }
    
          #toc {
            float: left;
            width: 100%;
          }
    
          #toc li {
            float: left;
            width: 33%;
          }
        </style>
      </head>
      <body>
        <div id=\"wrapper\">
          <h1>General Information</h1>
    
          <dl>
            <dt>Number of Pathways</dt>
            <dd>", length(kegg.collection), "</dd>
          </dl>
    
          <ul>
            <li>
              <a href=\"0verview Heatmaps.pdf\" target=\"_blank\">
                PSF signal Heatmaps (PDF)
              </a>
            </li>
          </ul>
    
       <ul id=\"toc\">", sep="", file=outfile)
    
    
      cat("
          </ul>
    
          <h1>Pathways</h1>
    
          <p>
            Population maps of the node and of the sink genes, and profile plots
            of the mean and maximum sink node signals are shown.
            Additionally the signal flow in each sample is given for the pathways.
          </p>", sep="", file=outfile)
    
    
        cat("
    
          <table>
            <thead>
              <tr>
                <th>Pathway</th>
                <th>Report sheet</th>
              </tr>
            </thead>
            <tbody>", sep="", file=outfile)
    
        for (ii in seq_along(kegg.collection))
        {
          cat("
              <tr>
                <td>", names(kegg.collection)[ii], "</td>
                <td><a href=\"", substring(make.names(names(kegg.collection)[ii]), 1, 100),
                  ".pdf\" target=\"_blank\">PDF</a></td>
              </tr>", sep="", file=outfile)
        }
    
        cat("
            </tbody>
          </table>", sep="", file=outfile)
      
    
      cat("
        </div>
      </body>
    </html>", sep="", file=outfile)
  
    close(outfile)
    
    
  }
  
}
