pipeline.htmlGenesetAnalysis <- function()
{
  if (!preferences$geneset.analysis)
  {
    return()
  }

  filename <- file.path(paste(files.name, "- Results"), "Geneset Analysis", "0verview.html")
  util.info("Writing:", filename)
  outfile <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Geneset Analysis Summary of ", files.name, " dataset</title>
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
    </style>
  </head>
  <body>
    <div id=\"wrapper\">
      <h1>General Information</h1>

      <dl>
        <dt>Number of Genesets</dt>
        <dd>", length(gs.def.list), "</dd>
        <dt>Categories</dt>
        <dd>", paste(paste(names(table(gs.def.list.categories)),
               table(gs.def.list.categories), sep=" ("), collapse=") , ") , ")</dd>
        <dt>Table of all GSZ scores</dt>
        <dd>
          <a href=\"../CSV Sheets/Sample GSZ scores.csv\" target=\"_blank\">
            Sample GSZ Scores.csv
          </a>
        </dd>
        <dt>GSZ/Fisher analysis</dt>
        <dd>
          <a href=\"0verview Heatmaps.pdf\" target=\"_blank\">
            Overview Heatmaps.pdf
          </a>
        </dd>
        <dt>Cancer Hallmarks</dt>
        <dd>
          <a href=\"0verview Cancer Hallmarks.pdf\" target=\"_blank\">
            Overview Cancer Hallmarks.pdf
          </a>
        </dd>
        <dt>Chromosome Expression</dt>
        <dd>
          <a href=\"0verview Chromosome Expression.pdf\" target=\"_blank\">
            Overview Chromosome Expression.pdf
          </a>
        </dd>
      </dl>

      <h2>Quick Links</h2>

      <ul>", sep="", file=outfile)

    for (i in names(table(gs.def.list.categories)))
    {
      cat("
        <li><a href=\"#", i, "\">", i, "</a></li>", sep="", file=outfile)
    }

  cat("
      </ul>

      <h1>Gene Sets</h1>

      <p>
        Enrichment profiles of individual predefined gene sets are shown as
        bar plots across all samples. Additionally the log FC-expression
        profiles of the leading metagenes are shown.
        Further, members of each gene set are shown as population maps and
        listed in Excel-files.
      </p>", sep="", file=outfile)

  for (i in names(table(gs.def.list.categories)))
  {
    category.gs.list <- gs.def.list[which(gs.def.list.categories == i)]

    cat("

      <h2 id=\"", i, "\">Category ", i, "</h2>

      <table>
        <thead>
          <tr>
            <th>Geneset name</th>
            <th>Category</th>
            <th>Profile</th>
            <th>Population Map</th>
            <th>Members</th>
          </tr>
        </thead>
        <tbody>", sep="", file=outfile)

    for (ii in 1:length(category.gs.list))
    {
      cat("
          <tr>
            <td>", names(category.gs.list)[ii], "</td>
            <td>", sapply(category.gs.list,function(x){x$Type})[ii], "</td>
            <td><a href=\"", substring(make.names(names(category.gs.list)[ii]), 1, 100),
              " profile.pdf\" target=\"_blank\">PDF</a></td>
            <td><a href=\"", substring(make.names(names(category.gs.list)[ii]), 1, 100),
              " map.pdf\" target=\"_blank\">PDF</a></td>
            <td><a href=\"", substring(make.names(names(category.gs.list)[ii]), 1, 100),
              ".csv\" >CSV</a></td>
          </tr>", sep="", file=outfile)
    }

    cat("
        </tbody>
      </table>", sep="", file=outfile)
  }

  cat("
    </div>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
