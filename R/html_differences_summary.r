pipeline.htmlDifferencesSummary <- function()
{
  dirname <- file.path(paste(files.name, "- Results"),
                       "Summary Sheets - Differences")

  if (!file.exists(dirname))
  {
    return()
  }

  filename <- file.path(dirname, "0verview.html")
  util.info("Writing:", filename)
  f <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Differences Summary of ", files.name, " dataset</title>
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
        display: block;
      }
    </style>
  </head>
  <body>
    <div id=\"wrapper\">
      <h1>Differences Analyses</h1>

      <h2>Pairwise Difference Summary Sheets</h2>

      <p>
        For each pairwise comparison a report sheet, gene lists and gene set
        lists are created analogous to the Sample Summary Sheets.
      </p>

      <table>
        <thead>
          <tr>
            <th>Comparison</th>
            <th>Summary Sheet</th>
            <th>Global Gene List</th>
            <th>Gene Set List</th>
          </tr>
        </thead>
        <tbody>", sep="", file=f)

  for (m in 1:ncol(indata))
  {
    name <- colnames(indata)[m]
    fname <- make.names(name)

    cat("
          <tr>
            <td>", name, "</td>
            <td>
              <a href=\"Reports/", fname, ".pdf\" target=\"_blank\">
                PDF
              </a>
            </td>
            <td>
              <a href=\"CSV Sheets/Gene Lists - Global/", fname, ".csv\" target=\"_blank\">
                CSV
              </a>
            </td>
            <td>
              <a href=\"CSV Sheets/Gene Set Lists - Global/", fname, ".csv\" target=\"_blank\">
                CSV
              </a>
            </td>
          </tr>", sep="", file=f)
  }

  cat("
        </tbody>
      </table>
    </div>
  </body>
</html>", sep="", file=f)

  close(f)
}
