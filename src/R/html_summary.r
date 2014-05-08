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
    </style>
  </head>
  <body>
    <div id=\"wrapper\">
      <h1>General Information</h1>

      <dl>
        <dt>Number of samples</dt>
        <dd>", ncol(indata), " (all replicates: ", ncol(indata.original), ")</dd>
        <dt>Number of categories</dt>
        <dd>", length(unique(group.labels)), "</dd>
        <dt>Number of genes</dt>
        <dd>", nrow(indata), "</dd>
        <dt>Dimension 1st level SOM</dt>
        <dd>", preferences$dim.som1, " x ", preferences$dim.som1, "</dd>
        <dt>Analysis finished</dt>
        <dd>", format(Sys.time(), "%a %b %d %X %Y %Z"), "</dd>
      </dl>

      <h1>Results</h1>

      <h2>1st Level SOM Analysis</h2>

      <p>
        These reports show the collection of first level SOM images of all
        tissue samples. Alternative SOM images are shown with different color
        scales and also as rank SOM using FC-, WAD-, and Shrinkage-t-scores.
        Further, supporting maps, entropy and topology profiles provide
        supplementary information about the 1st level SOM.
      </p>

      <ul>
        <li>
          <a href=\"Expression Portraits.pdf\" target=\"_blank\">
            1st Level SOM Expression Images
          </a>
        </li>
        <li>
          <a href=\"Expression Portraits alternative.pdf\" target=\"_blank\">
            Alternative Color Scales
          </a>
        </li>
        <li>
          <a href=\"Rank Maps.pdf\" target=\"_blank\">
            Rank Profiles
          </a>
        </li>
      </ul>

      <ul>
        <li>
          <a href=\"Summary Sheets - Groups/Expression Portraits Groups.pdf\" target=\"_blank\">
            1st Level SOM Expression Images - Group Specific
          </a>
        </li>
        <li>
          <a href=\"Supporting Maps&amp;Profiles/Supporting Maps.pdf\" target=\"_blank\">
            Supporting Maps
          </a>
        </li>
      </ul>

      <ul>
        <li>
          <a href=\"Supporting Maps&amp;Profiles/Entropy Profiles.pdf\" target=\"_blank\">
            Entropy Profiles
          </a>
        </li>
        <li>
          <a href=\"Supporting Maps&amp;Profiles/Topology Profiles.pdf\" target=\"_blank\">
            Topology Profiles
          </a>
        </li>
      </ul>

      <h2>Sample Summaries</h2>

      <p>
        For each sample a report sheet is created which summarizes the most
        relevant information.
      </p>

      <ul>
        <li>
          <a href=\"Summary Sheets - Samples/0verview.html\" target=\"_blank\">
            Sample Report Sheets
          </a>
        </li>
      </ul>", sep="", file=outfile)

  if (preferences$geneset.analysis)
  {
    cat("
      <h2>Geneset Enrichment Analysis</h2>

      <p>
        Enrichment of predefined gene sets in the samples is calculated and
        visualized as heatmaps, profile plots and population maps.
      </p>

      <ul>
        <li>
          <a href=\"Geneset Analysis/0verview.html\" target=\"_blank\">
            Geneset Analysis Results
          </a>
        </li>
      </ul>", sep="", file=outfile)
  }

  cat("
      <h2>2nd Level Metagene Analysis</h2>

      <p>
        Several agglomerative methods based either on distance or on
        correlation metrics are applied to the samples using filtered subsets
        of metagenes.
        The reports show cluster heatmaps and dendrograms, pairwise
        correlation maps, correlation spanning trees and ICA results.
      </p>

      <ul>
        <li>
          <a href=\"2nd lvl Metagene Analysis/2nd lvl SOM.pdf\" target=\"_blank\">
            2nd Level SOM
          </a>
        </li>
        <li>
          <a href=\"2nd lvl Metagene Analysis/Similarity Analysis.pdf\" target=\"_blank\">
            Similarity Based Methods (Neighbor Joining &amp; Hierarchical Clustering)
          </a>
        </li>
        <li>
          <a href=\"2nd lvl Metagene Analysis/Correlation Analysis.pdf\" target=\"_blank\">
            Correlation Based Methods (MST, Correlation Graphs, PCM)
          </a>
        </li>
        <li>
          <a href=\"2nd lvl Metagene Analysis/Component Analysis.pdf\" target=\"_blank\">
            Component Based Methods (2d-ICA, 3d-ICA)
          </a>
        </li>
      </ul>

      <h2>3rd Level Spot Analysis</h2>

      <p>
        Different criteria of spot selection such as overexpression or
        mutual correlations between the metagenes where applied.
        Summary sheets of integral SOM maps and subsequent analyses base on
        this third level information.
      </p>

      <ul>
        <li>
          <a href=\"Summary Sheets - Integral/0verview.html\" target=\"_blank\">
            Spot Report Sheets
          </a>
        </li>
      </ul>
    </div>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
