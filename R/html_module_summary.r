pipeline.htmlModuleSummary <- function(env)
{
  filename <- file.path("Summary Sheets - Modules","0verview.html")

  util.info("Writing:", filename)
  outfile <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Spot Summary of ", env$files.name, " dataset</title>
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

      table a {
        margin: 0 6px;
      }
    </style>
  </head>
  <body>
    <div id=\"wrapper\">", sep="", file=outfile)


  #### Report Sheets ####

  cat("<h1>Expression Module Reports</h1>

      <p>
        Expression module detection applies different criteria to define spot clusters of co-regulated metagenes:
        overexpression and underexpression, mutual correlations between the metagenes,
        k-means clustering and cluster detection in the distance map.<br>

        <ul>
          <li>
          General report sheets provide overview maps of the modules detected, module-related heatmaps characterizing the
          expression profiles of the selected features in the series of samples, and individual module reports
          providing detailed information such as ranked gene list selected with the module and functional analysis in terms of
          gene set lists ranked for over-representation in the module.
          </li>
        </ul>
        <ul>
          <li>
          Profile sheets provide the module expression profiles on single sample and group resolution.
          </li>
        </ul>
        <ul>
          <li>
          Module relation reports contain the wTO correlation network and association of the groups to the modules.
          </li>
        </ul>
        <ul>
          <li>
          Chromosome plots visualize enrichment of chromosomal postitions (chromosome/band) in the modules as overview heatmaps and individual chromosome plots. 
          </li>
        </ul>
      </p>", sep="", file=outfile)

  
  cat("<h2>Overexpression Spots</h2>

      <ul>
        <li>
          <a href=\"Overexpression Spots/Report.pdf\" target=\"_blank\">
                General module report (PDF)
              </a>
        </li>
        <li>
          <a href=\"Overexpression Spots/Profiles.pdf\" target=\"_blank\">
                Module Profiles (PDF)
              </a>
        </li>
        <li>
          <a href=\"Overexpression Spots/Relations.pdf\" target=\"_blank\">
                Module Relation (PDF)
              </a>
        </li>
        <li>
          <a href=\"Overexpression Spots/Chromosomes.pdf\" target=\"_blank\">
                Module Enrichment in Chromosomes (PDF)
              </a>
        </li>
      </ul>", sep="", file=outfile)

  if (length(unique(env$group.labels)) > 1)
  {
    cat("<h2>Group Overexpression Spots</h2>
  
        <ul>
          <li>
            <a href=\"Group Overexpression Spots/Report.pdf\" target=\"_blank\">
                  General module report (PDF)
                </a>
          </li>
          <li>
            <a href=\"Group Overexpression Spots/Profiles.pdf\" target=\"_blank\">
                  Module Profiles (PDF)
                </a>
          </li>
          <li>
            <a href=\"Group Overexpression Spots/Relations.pdf\" target=\"_blank\">
                  Module Relation (PDF)
                </a>
          </li>
          <li>
            <a href=\"Group Overexpression Spots/Chromosomes.pdf\" target=\"_blank\">
                  Module Enrichment in Chromosomes (PDF)
                </a>
          </li>
        </ul>", sep="", file=outfile)
  }

  cat("<h2>Underexpression Spots</h2>

      <ul>
        <li>
          <a href=\"Underexpression Spots/Report.pdf\" target=\"_blank\">
                General module report (PDF)
              </a>
        </li>

      </ul>", sep="", file=outfile)

  cat("<h2>Correlation Cluster</h2>

      <ul>
        <li>
          <a href=\"Correlation Cluster/Report.pdf\" target=\"_blank\">
                General module report (PDF)
              </a>
        </li>

      </ul>", sep="", file=outfile)
  
  cat("<h2>K-Means Cluster</h2>

      <ul>
        <li>
          <a href=\"K-Means Cluster/Report.pdf\" target=\"_blank\">
                General module report (PDF)
              </a>
        </li>
        <li>
          <a href=\"K-Means Cluster/Profiles.pdf\" target=\"_blank\">
                Module Profiles (PDF)
              </a>
        </li>
        <li>
          <a href=\"K-Means Cluster/Relations.pdf\" target=\"_blank\">
                Module Relation (PDF)
              </a>
        </li>
        <li>
          <a href=\"K-Means Cluster/Chromosomes.pdf\" target=\"_blank\">
                Module Enrichment in Chromosomes (PDF)
              </a>
        </li>
      </ul>", sep="", file=outfile)
  
  cat("<h2>D-Cluster</h2>

      <ul>
        <li>
          <a href=\"D-Cluster/Report.pdf\" target=\"_blank\">
                General module report (PDF)
              </a>
        </li>
        <li>
          <a href=\"D-Cluster/Profiles.pdf\" target=\"_blank\">
                Module Profiles (PDF)
              </a>
        </li>
        <li>
          <a href=\"D-Cluster/Relations.pdf\" target=\"_blank\">
                Module Relation (PDF)
              </a>
        </li>
        <li>
          <a href=\"D-Cluster/Chromosomes.pdf\" target=\"_blank\">
                Module Enrichment in Chromosomes (PDF)
              </a>
        </li>
      </ul>", sep="", file=outfile)  
  
  #### CSV Sheets ####

  cat("<h1>Expression Module CSV Sheets</h1>

      <table>
        <thead>
          <tr>
            <th></th>
            <th>CSV Tables</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>
              Overexpression Spots
            </td>
            <td>", sep="", file=outfile)

  for (m in seq_along(env$spot.list.overexpression$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/Overexpression Spots ",
        names(env$spot.list.overexpression$spots)[m],".csv\" target=\"_blank\">",
        names(env$spot.list.overexpression$spots)[m],"</a>", sep="", file=outfile)
  }

  if (length(unique(env$group.labels)) > 1)
  {
    cat("</td>
          </tr>
          <tr>
            <td>
                Group Overexpression Spots
            </td>
            <td>", sep="", file=outfile)

    for (m in seq_along(env$spot.list.group.overexpression$spots))
    {
      cat("<a href=\"../CSV Sheets/Spot Lists/Group Overexpression Spots ",
          names(env$spot.list.group.overexpression$spots)[m],".csv\" target=\"_blank\">",
          names(env$spot.list.group.overexpression$spots)[m],"</a>", sep="", file=outfile)
    }

  }

  cat("</td>
          </tr>
          <tr>
            <td>
                K-Means Cluster
            </td>
            <td>", sep="", file = outfile)

  for (m in seq_along(env$spot.list.kmeans$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/K-Means Cluster ",
        names(env$spot.list.kmeans$spots)[m],".csv\" target=\"_blank\">",
        names(env$spot.list.kmeans$spots)[m],"</a>", sep="", file=outfile)
  }

  
  cat("</td>
          </tr>
          <tr>
            <td>
                Distance Map Cluster
            </td>
            <td>", sep="", file = outfile)
  
  for (m in seq_along(env$spot.list.dmap$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/D-Cluster ",
        names(env$spot.list.dmap$spots)[m],".csv\" target=\"_blank\">",
        names(env$spot.list.dmap$spots)[m],"</a>", sep="", file=outfile)
  }
  
  cat("</td>
          </tr>
        </tbody>
      </table>
    </div>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
