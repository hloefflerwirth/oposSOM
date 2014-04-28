pipeline.htmlSummary <- function()
{
  filename <- file.path(paste(files.name, "- Results"), "Summary.html")
  util.info("Writing:", filename)
  outfile <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Summary of ", files.name, " dataset</title>

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
    </style>
  </head>
  <body>
    <h1>General Information</h1>

    <table>
      <tr>
        <td>Number of samples:</td>
        <td>", ncol(indata), " (all replicates: ", ncol(indata.original), ")</td>
      </tr>
      <tr>
        <td>Number of categories:</td>
        <td>", length(unique(group.labels)), "</td>
      </tr>
      <tr>
        <td>Number of genes:</td>
        <td>", nrow(indata), "</td>
      </tr>
      <tr>
        <td>Analysis finished:</td>
        <td>", format(Sys.time(), "%a %b %d %X %Y %Z"), "</td>
      </tr>
    </table>

    <table>
      <tr>
        <td>Dimension 1st level SOM:</td>
        <td>", preferences$dim.som1, "x", preferences$dim.som1, "</td>
      </tr>
      <tr>
        <td>Dimension 2nd level SOM:</td>
        <td>", preferences$dim.som2, "x", preferences$dim.som2, "</td>
      </tr>
    </table>

    <h1>Results</h1>

    <table>
      <thead>
        <tr>
          <th colspan=\"2\">1st Level SOM Analysis</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td rowspan=\"5\">
            These reports show the collection of first level SOM images of all
            tissue samples. Alternative SOM images are shown with different color
            scales and also as rank SOM using FC-, WAD-, and Shrinkage-t-scores.
            Further, supporting maps, entropy and topology profiles provide
            supplementary information about the 1st level SOM.
          </td>
          <td>
            <a href=\"Expression Portraits.pdf\" target=\"_blank\">
              1st level SOM expression images
            </a>
            <a href=\"Expression Portraits alternative.pdf\" target=\"_blank\">
              alternative color scales
            </a>
            <a href=\"Rank Maps.pdf\" target=\"_blank\">
              rank profiles
            </a>
          </td>
        </tr>
        <tr>
          <td>
            <a href=\"Summary Sheets - Groups/Expression Portraits Groups.pdf\" target=\"_blank\">
              1st level SOM expression images - group specific
            </a>
          </td>
        </tr>
        <tr>
          <td>
            <a href=\"Supporting Maps&Profiles/Supporting Maps.pdf\" target=\"_blank\">
              Supporting Maps
            </a>
          </td>
        </tr>
        <tr>
          <td>
            <a href=\"Supporting Maps&Profiles/Entropy Profiles.pdf\" target=\"_blank\">
              Entropy Profiles
            </a>
            <a href=\"Supporting Maps&Profiles/Topology Profiles.pdf\" target=\"_blank\">
              Topology Profiles
            </a>
          </td>
        </tr>
      </tbody>
    </table>

    <table>
      <thead>
        <tr>
          <th colspan=\"5\">Sample Summaries</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td rowspan=\"5\">
            For each sample a report sheet is created which summarizes the most
            relevant information.
          </td>
          <td>
            <a href=\"Summary Sheets - Samples/0verview.html\" target=\"_blank\">
              browse sample report sheets
            </a>
          </td>
        </tr>
      </tbody>
    </table>", sep="", file=outfile)

  if (preferences$geneset.analysis)
  {
    cat("
    <table>
      <thead>
        <tr>
          <th colspan=\"5\">
            Gene Set Enrichment Analysis
          </th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td rowspan=\"5\">
            Enrichment of predefined gene sets in the samples is calculated and
            visualized as heatmaps, profile plots and population maps.
          </td>
          <td>
            <a href=\"Geneset Analysis/0verview.html\" target=\"_blank\">
              browse geneset analysis results
            </a>
          </td>
        <tr>
      </tbody>
    </table>", sep="", file=outfile)
  }

  cat("
    <table>
      <thead>
        <tr>
          <th colspan=\"4\">2nd Level Metagene Analysis</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td rowspan=\"5\">
            Several agglomerative methods based either on distance or on
            correlation metrics are applied to the samples using filtered subsets
            of metagenes.
            The reports show cluster heatmaps and dendrograms, pairwise
            correlation maps, correlation spanning trees and ICA results.
          </td>
          <td>
            <a href=\"2nd lvl Metagene Analysis/2nd lvl SOM.pdf\" target=\"_blank\">
              Second level SOM
            </a>
          </td>
        </tr>
        <tr>
          <td>
            <a href=\"2nd lvl Metagene Analysis/Similarity Analysis.pdf\" target=\"_blank\">
              Similarity based methods (Neighbor Joining & Hierarchical Clustering)
            </a>
          </td>
        </tr>
        <tr>
          <td>
            <a href=\"2nd lvl Metagene Analysis/Correlation Analysis.pdf\" target=\"_blank\">
              Correlation based methods (MST, Correlation Graphs, PCM)
            </a>
          </td>
        </tr>
        <tr>
          <td>
            <a href=\"2nd lvl Metagene Analysis/Component Analysis.pdf\" target=\"_blank\">
              Component based methods (2d-ICA, 3d-ICA)
            </a>
          </td>
        </tr>
      </tbody>
    </table>

    <table>
      <thead>
        <tr>
          <th colspan=\"5\">3rd Level Spot Analysis</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td rowspan=\"6\">
            Different criteria of spot selection such as overexpression or
            mutual correlations between the metagenes where applied.
            Summary sheets of integral SOM maps and subsequent analyses base on
            this third level information.
          </td>
          <td>
            <a href=\"Summary Sheets - Integral/0verview.html\" target=\"_blank\">
              browse spot report sheets
            </a>
          </td>
        </tr>
      </tbody>
    </table>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
