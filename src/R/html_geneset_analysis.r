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
    <title>Geneset Analysis Summary of ", files.name, "dataset</title>
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
    <h1>General Information</h1>

    <table>
      <tr>
        <td>Number of Genesets:</td>
        <td>", length(gs.def.list), "</td>
      </tr>
      <tr>
        <td>Categories:</td>
        <td>", paste(paste(names(table(gs.def.list.categories)),
               table(gs.def.list.categories), sep=" ("), collapse=") , ") , ")</td>
      </tr>
      <tr>
        <td>Table of all GSZ scores:</td>
        <td><a href=\"../CSV sheets/Sample GSZ scores.csv\" target=\"_blank\">CSV</a></td>
      </tr>
      <tr>
        <td>GSZ/Fisher analysis: Heatmaps & p-value histograms:</td>
        <td><a href=\"0verview Heatmaps.pdf\" target=\"_blank\">PDF</a></td>
      </tr>
    </table>

    <h1>Quick Links</h1>

    <p>Go to category:", sep="", file=outfile)

  for (i in names(table(gs.def.list.categories)))
  {
    cat("<a href=\"#", i, "\">", i, "</a>", sep="", file=outfile)
  }

  cat("</p>

    <h1>Gene Sets</h1>

    <table>
      <tr>
        <td class=\"justy\">
          Enrichment profiles of individual predefined gene sets are shown as
          bar plots across all samples. Additionally the log FC-expression
          profiles of the leading metagenes are shown.
          Further, members of each gene set are shown as population maps and
          listed in Excel-files.
        </td>
      </tr>
    </table>", sep="", file=outfile)

  for (i in names(table(gs.def.list.categories)))
  {
    category.gs.list <- gs.def.list[which(gs.def.list.categories == i)]

    cat("

    <h2 id=\"", i, "\">Category ", i, "</h2>

    <table>
      <thead>
        <tr>
          <th style=\"width: 55%\">Geneset name</th>
          <th style=\"width: 9%\">Category</th>
          <th style=\"width: 12%\">Profile</th>
          <th style=\"width: 12%\">Population Map</th>
          <th style=\"width: 12%\">Members</th>
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
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
