pipeline.htmlIntegralSummary <- function()
{
  filename <- file.path(paste(files.name, "- Results"), "Summary Sheets - Integral", "0verview.html")
  util.info("Writing:", filename)
  outfile <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Spot Summary of ", files.name, "dataset</title>
    <style>
      body {
        background-color: #fff;
      }

      p {
        margin-top: 1px;
        margin-bottom: 1px;
        padding-left: 10px;
        text-indent: -10px
      }

      table {
        width: 90%;
        margin: 40px auto;
        border: 2px solid #333;
      }

      table td {
        width: 50%;
      }

      table a {
        margin: 0 10px;
      }

      .justy {
        text-align: justify;
      }
    </style>
  </head>
  <body>
    <h1>Integral Maps Summary Sheets</h1>

    <table>
      <tr>
        <td class=\"justy\">
          These analyses apply different criteria of spot selection such as
          overexpression, underexpression, maximum and minimum of metagene
          expression. Additionaly, mutual correlations between the metagenes
          and k-means clustering was applied to define spots of co-regulated
          metagenes.
          Gene set enrichment analysis provides the leading sets of each of
          the spots considered. Spot-related heatmaps characterize the
          expression profiles of the selected features in the series of samples.
          Single spot summary sheets provide detailed information on each of
          the spots such as the ranked list of samples which overexpress this
          feature in decreasing order according to the mean t-shrinkage
          statistics of the spot and the ranked list of the overrepresented
          gene sets together with the histogram of the respective p-value
          distribution.
        </td>
      </tr>
    </table>

    <table border=2, width=90%>
      <thead>
        <tr>
          <th style=\"width: 25%\">Type</th>
          <th style=\"width: 20%\">Summary Sheets</th>
          <th style=\"width: 55%\">Spot's CSV tables</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>Overexpression</td>
          <td><a href=\"Overexpression.pdf\" target=\"_blank\">PDF</a></td>
          <td>", sep="", file=outfile)

  for (m in 1:length(GS.infos.overexpression$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/Sample-Overexpression ",
        names(GS.infos.overexpression$spots)[m],".csv\" target=\"_blank\">",
        names(GS.infos.overexpression$spots)[m],"</a>", sep="", file=outfile)
  }

  cat("</td>
        </tr>
        <tr>
          <td>Underexpression</td>
          <td><a href=\"Underexpression.pdf\" target=\"_blank\">PDF</a></td>
          <td>", sep="", file=outfile)

  for (m in 1:length(GS.infos.underexpression$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/Sample-Underexpression ",
        names(GS.infos.underexpression$spots)[m],".csv\" target=\"_blank\">",
        names(GS.infos.underexpression$spots)[m], "</a>", sep="", file=outfile)
  }

  cat("</td>
        </tr>
        <tr>
          <td>Mutual Correlation</td>
          <td><a href=\"Correlation Cluster.pdf\" target=\"_blank\">PDF</a></td>
          <td>", sep="", file=outfile)

  for (m in 1:length(GS.infos.correlation$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/Correlation Cluster ",
        names(GS.infos.correlation$spots)[m],".csv\" target=\"_blank\">",
        names(GS.infos.correlation$spots)[m],"</a>", sep="", file=outfile)
  }

  cat("</td>
        </tr>
        <tr>
          <td>k-Means Clustering</td>
          <td><a href=\"K-Means Cluster.pdf\" target=\"_blank\">PDF</a></td>
          <td>", sep="", file = outfile)

  for (m in 1:length(GS.infos.kmeans$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/K-Means Cluster ",
        names(GS.infos.kmeans$spots)[m],".csv\" target=\"_blank\">",
        names(GS.infos.kmeans$spots)[m],"</a>", sep="", file=outfile)
  }

  cat("</td>
        </tr>
        <tr>
          <td>Group Overexpression</td>
          <td><a href=\"Group Overexpression.pdf\" target=\"_blank\">PDF</a></td>
          <td>", sep="", file=outfile)

  for (m in 1:length(GS.infos.group.overexpression$spots))
  {
    cat("<a href=\"../CSV Sheets/Spot Lists/Group Overexpression ",
        names(GS.infos.group.overexpression$spots)[m],".csv\" target=\"_blank\">",
        names(GS.infos.group.overexpression$spots)[m],"</a>", sep="", file=outfile)
  }

  cat("</td>
        </tr>
      </tbody>
    </table>

    <h1>Spot Report Sheets</h1>

    <table>
      <tr>
        <td class=\"justy\">
          Spot repots contain assignment of the spots to featuring samples and
          to respective sample groups.
          Additionally, relative numbers of significant gene sets are given for
          each spot and gene set category.
        </td>
      </tr>
    </table>

    <table>
      <thead>
        <tr>
          <th>Type</th>
          <th>Summary Sheets</th>
        </tr>
      </thead>
      <tbody>
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
      </tbody>
    </table>

    <h1>Spot Network Analysis</h1>

    <table>
      <tr>
        <td class=\"justy\">
          Networks of spot association are visualized as graphs. Both,
          correlation networks (connecting spots with r>0.5) and correlation
          spanning trees, are given for individual spots and also for patterns
          of spots occuring in the sample SOM images.
        </td>
      </tr>
    </table>

    <table>
      <thead>
        <tr>
          <th>Type</th>
          <th>Summary Sheets</th>
        </tr>
      </thead>
      <tbody>
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
      </tbody>
    </table>

    <h1>Chromosomal Enrichment</h1>

    <table>
      <tr>
        <td class=\"justy\">
          For each spot, enrichment of chromosomal postitions (chromosome/band)
          is visualized as overview heatmaps and individual chromosome plots.
        </td>
      </tr>
    </table>

    <table>
      <thead>
        <tr>
          <th>Type</th>
          <th>Summary Sheets</th>
        </tr>
      </thead>
      <tbody>
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
      </tbody>
    </table>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
