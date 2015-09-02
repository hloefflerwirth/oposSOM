pipeline.htmlIntegralSummary <- function()
{
  filename <- file.path(paste(files.name, "- Results"),
                        "Summary Sheets - Integral",
                        "0verview.html")

  util.info("Writing:", filename)
  outfile <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Spot Summary of ", files.name, " dataset</title>
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
    <div id=\"wrapper\">
      <h1>Integral Maps Summary Sheets</h1>

      <p>
        These analyses apply different criteria to define spot clusters of co-regulated metagenes:
        overexpression and underexpression, mutual correlations between the metagenes,
        k-means clustering and cluster detection in the distance map.
        Gene set enrichment analysis provides the leading sets of each of
        the spots considered. Spot-related heatmaps characterize the
        expression profiles of the selected features in the series of samples.
        Single spot summary sheets provide detailed information on each of
        the spots such as the ranked list of samples which overexpress this
        feature in decreasing order according to the mean t-statistics of the
        spot and the ranked list of the overrepresented
        gene sets together with the histogram of the respective p-value
        distribution.
      <p>

      <table>
        <thead>
          <tr>
            <th>Summary Sheets</th>
            <th>Spot's CSV tables</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>
              <a href=\"Overexpression.pdf\" target=\"_blank\">
                Overexpression (PDF)
              </a>
            </td>
            <td>", sep="", file=outfile)

  for (m in seq_along(spot.list.overexpression$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/Overexpression Spots ",
        names(spot.list.overexpression$spots)[m],".csv\" target=\"_blank\">",
        names(spot.list.overexpression$spots)[m],"</a>", sep="", file=outfile)
  }

  cat("</td>
          </tr>
          <tr>
            <td>
              <a href=\"../3rd lvl Spot Analysis/Signature Sets.csv\" target=\"_blank\">
                Overexpresion Signatures (CSV)
              </a>
            </td>
            <td></td>
          </tr>", sep="", file=outfile)

  if (length(unique(group.labels)) > 1)
  {
    cat("
          <tr>
            <td>
              <a href=\"Group Overexpression.pdf\" target=\"_blank\">
                Group Overexpression (PDF)
              </a>
            </td>
            <td>", sep="", file=outfile)

    for (m in seq_along(spot.list.group.overexpression$spots))
    {
      cat("<a href=\"../CSV Sheets/Spot Lists/Group Overexpression Spots ",
          names(spot.list.group.overexpression$spots)[m],".csv\" target=\"_blank\">",
          names(spot.list.group.overexpression$spots)[m],"</a>", sep="", file=outfile)
    }

    cat("</td>
          </tr>", sep="", file=outfile)
  }

  cat("
          <tr>
            <td>
              <a href=\"Underexpression.pdf\" target=\"_blank\">
                Underexpression (PDF)
              </a>
            </td>
            <td>", sep="", file=outfile)

  for (m in seq_along(spot.list.underexpression$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/Underexpression Spots ",
        names(spot.list.underexpression$spots)[m],".csv\" target=\"_blank\">",
        names(spot.list.underexpression$spots)[m], "</a>", sep="", file=outfile)
  }

  cat("</td>
          </tr>
          <tr>
            <td>
              <a href=\"Correlation Cluster.pdf\" target=\"_blank\">
                Correlation Cluster (PDF)
              </a>
            </td>
            <td>", sep="", file=outfile)

  for (m in seq_along(spot.list.correlation$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/Correlation Clusters ",
        names(spot.list.correlation$spots)[m],".csv\" target=\"_blank\">",
        names(spot.list.correlation$spots)[m],"</a>", sep="", file=outfile)
  }

  cat("</td>
          </tr>
          <tr>
            <td>
              <a href=\"K-Means Cluster.pdf\" target=\"_blank\">
                K-Means Cluster (PDF)
              </a>
            </td>
            <td>", sep="", file = outfile)

  for (m in seq_along(spot.list.kmeans$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/K-Means Clusters ",
        names(spot.list.kmeans$spots)[m],".csv\" target=\"_blank\">",
        names(spot.list.kmeans$spots)[m],"</a>", sep="", file=outfile)
  }

  
  cat("</td>
          </tr>
          <tr>
            <td>
              <a href=\"D-Clusters.pdf\" target=\"_blank\">
                Distance Map Clusters (PDF)
              </a>
            </td>
            <td>", sep="", file = outfile)
  
  for (m in seq_along(spot.list.dmap$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/D-Clusters ",
        names(spot.list.dmap$spots)[m],".csv\" target=\"_blank\">",
        names(spot.list.dmap$spots)[m],"</a>", sep="", file=outfile)
  }
  
  cat("</td>
          </tr>
        </tbody>
      </table>

      <h1>Spot Module Report Sheets</h1>

      <p>
        Reports contain the spot module expression profiles and assignments of
        the spots to samples and to groups.
      </p>

      <ul>
        <li>
          <a href=\"../3rd lvl Spot Analysis/Spot Report - Overexpression Spots.pdf\" target=\"_blank\">
            Overexpression Spot Report (PDF)
          </a>
        </li>
        <li>
          <a href=\"../3rd lvl Spot Analysis/Spot Report - Underexpression Spots.pdf\" target=\"_blank\">
            Underexpression Spot Report (PDF)
          </a>
        </li>
        <li>
          <a href=\"../3rd lvl Spot Analysis/Spot Report - K-Means Clusters.pdf\" target=\"_blank\">
            K-Means Cluster Report (PDF)
          </a>
        </li>
        <li>
          <a href=\"../3rd lvl Spot Analysis/Spot Report - D-Clusters.pdf\" target=\"_blank\">
            D-Cluster Report (PDF)
          </a>
        </li>", sep="", file=outfile)

  if (length(unique(group.labels)) > 1)
  {
    cat("
        <li>
          <a href=\"../3rd lvl Spot Analysis/Spot Report - Group Overexpression Spots.pdf\" target=\"_blank\">
            Group Overexpression Report (PDF)
          </a>
        </li>", sep="", file=outfile)
  }

  cat("
      </ul>

      <h1>Spot Module Network Analysis</h1>

      <p>
        Networks of spot association are visualized as graphs. WTO,
        correlation networks and correlation spanning trees, are given for
        individual spots and spot patterns.
      </p>

      <ul>
        <li>
          <a href=\"../3rd lvl Spot Analysis/wTO Networks - Overexpression Spots.pdf\" target=\"_blank\">
            Overexpression Networks (PDF)
          </a>
        </li>
        <li>
          <a href=\"../3rd lvl Spot Analysis/wTO Networks - Underexpression Spots.pdf\" target=\"_blank\">
            Underexpression Networks (PDF)
          </a>
        </li>
        <li>
          <a href=\"../3rd lvl Spot Analysis/wTO Networks - K-Means Clusters.pdf\" target=\"_blank\">
            K-Means Cluster Networks (PDF)
          </a>
        </li>
        <li>
          <a href=\"../3rd lvl Spot Analysis/wTO Networks - D-Clusters.pdf\" target=\"_blank\">
            D-Cluster Networks (PDF)
          </a>
        </li>", sep="", file=outfile)

  if (length(unique(group.labels)) > 1)
  {
    cat("
        <li>
          <a href=\"../3rd lvl Spot Analysis/wTO Networks - Group Overexpression Spots.pdf\" target=\"_blank\">
            Group Overexpression Networks (PDF)
          </a>
        </li>", sep="", file=outfile)
  }

  cat("
      </ul>

      <h1>Chromosomal Enrichment</h1>

      <p>
        For each spot, enrichment of chromosomal postitions (chromosome/band)
        is visualized as overview heatmaps and individual chromosome plots.
      </p>

      <ul>
        <li>
          <a href=\"../3rd lvl Spot Analysis/Chromosomal Enrichment - Overexpression Spots.pdf\" target=\"_blank\">
            Chromosomal Enrichment of Overexpression Spots (PDF)
          </a>
        </li>
        <li>
          <a href=\"../CSV Sheets/Chromosomal Enrichment/Overexpression Spots.csv\" target=\"_blank\">
            Chromosomal Enrichment of Overexpression Spots (CSV)
          </a>
        </li>", sep="", file=outfile)

  if (length(unique(group.labels)) > 1)
  {
    cat("
        <li>
          <a href=\"../3rd lvl Spot Analysis/Chromosomal Enrichment - Group Overexpression Spots.pdf\" target=\"_blank\">
            Chromosomal Enrichment of Group Overexpression Spots (PDF)
          </a>
        </li>
        <li>
          <a href=\"../CSV Sheets/Chromosomal Enrichment/Group Overexpression Spots.csv\" target=\"_blank\">
            Chromosomal Enrichment of Group Overexpression Spots (CSV)
          </a>
        </li>", sep="", file=outfile)
  }

  cat("

        <li>
          <a href=\"../3rd lvl Spot Analysis/Chromosomal Enrichment - Underexpression Spots.pdf\" target=\"_blank\">
            Chromosomal Enrichment of Underexpression Chromosomal Enrichment (PDF)
          </a>
        </li>
        <li>
          <a href=\"../CSV Sheets/Chromosomal Enrichment/Underexpression Spots.csv\" target=\"_blank\">
            Chromosomal Enrichment of Underexpression Spot Chromosome Map (CSV)
          </a>
        </li>
        <li>
          <a href=\"../3rd lvl Spot Analysis/Chromosomal Enrichment - K-Means Clusters.pdf\" target=\"_blank\">
            Chromosomal Enrichment of K-Means Cluster Chromosomal Enrichment (PDF)
          </a>
        </li>
        <li>
          <a href=\"../CSV Sheets/Chromosomal Enrichment/K-Means Clusters.csv\" target=\"_blank\">
            Chromosomal Enrichment of K-Means Cluster Spot Chromosome Map (CSV)
          </a>
        </li>
        <li>
          <a href=\"../3rd lvl Spot Analysis/Chromosomal Enrichment - D-Clusters.pdf\" target=\"_blank\">
            Chromosomal Enrichment of D-Cluster Chromosomal Enrichment (PDF)
          </a>
        </li>
        <li>
          <a href=\"../CSV Sheets/Chromosomal Enrichment/D-Clusters.csv\" target=\"_blank\">
            Chromosomal Enrichment of D-Cluster Spot Chromosome Map (CSV)
          </a>
        </li>
      </ul>
    </div>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
