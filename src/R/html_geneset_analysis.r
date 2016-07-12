pipeline.htmlGenesetAnalysis <- function()
{
  if (!preferences$activated.modules$geneset.analysis)
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
        <dt>Number of Gene Sets</dt>
        <dd>", length(gs.def.list), "</dd>
        <dt>Categories</dt>
        <dd>", paste(paste(names(table(sapply(gs.def.list, function(x) { x$Type }))),
               table(sapply(gs.def.list, function(x) { x$Type })), sep=" ("), collapse=") , ") , ")</dd>
      </dl>

      <ul>
        <li>
          <a href=\"../CSV Sheets/Sample GSZ scores.csv\" target=\"_blank\">
            Table of sample GSZ scores (CSV)
          </a>
        </li>
        <li>
          <a href=\"0verview Heatmaps.pdf\" target=\"_blank\">
            Enrichment Score Heatmaps (PDF)
          </a>
        </li>
        <li>
          <a href=\"0verview Cancer Hallmarks.pdf\" target=\"_blank\">
            Hallmarks of Cancer (PDF)
          </a>
        </li>
      </ul>

      <h2>Category Links</h2>

      <ul id=\"toc\">", sep="", file=outfile)

    for (i in names(table(sapply(gs.def.list, function(x) { x$Type }))))
    {
      cat("
        <li><a href=\"#", i, "\">", i, "</a></li>", sep="", file=outfile)
    }

  cat("
      </ul>

      <h1>Gene Sets</h1>

      <p>
        The report sheets show the enrichment (GSZ-) profile as
        bar plot across all samples, the heatmap of most variant gene set members
        and the mapping of all member genesinto the SOM space.
        Further, all members of each individual gene set are given in the CSV spreadsheets.
        tables.
      </p>", sep="", file=outfile)

  for (i in names(table(sapply(gs.def.list, function(x) { x$Type }))))
  {
    category.gs.list <- gs.def.list[which(sapply(gs.def.list, function(x) { x$Type }) == i)]

    cat("

      <h2 id=\"", i, "\">Category ", i, "</h2>

      <table>
        <thead>
          <tr>
            <th>Geneset name</th>
            <th>Category</th>
            <th>Report Sheet</th>
            <th>Members</th>
          </tr>
        </thead>
        <tbody>", sep="", file=outfile)

    for (ii in seq_along(category.gs.list))
    {
      cat("
          <tr>
            <td>", names(category.gs.list)[ii], "</td>
            <td>", sapply(category.gs.list,function(x){x$Type})[ii], "</td>
            <td><a href=\"", substring(make.names(names(category.gs.list)[ii]), 1, 100),
              ".pdf\" target=\"_blank\">PDF</a></td>
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
