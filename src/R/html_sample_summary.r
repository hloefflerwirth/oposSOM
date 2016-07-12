pipeline.htmlSampleSummary <- function()
{
  if(ncol(indata) >= 1000) return()
	
  filename <- file.path(paste(files.name, "- Results"),
                        "Summary Sheets - Samples",
                        "0verview.html")

  util.info("Writing:", filename)
  outfile <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Sample Summary of ", files.name, " dataset</title>
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

      table td a {
        display: block;
      }
    </style>
  </head>
  <body>
    <div id=\"wrapper\">
      <h1>Group Overview</h1>

      <table>
        <thead>
          <tr>
            <th>Groups</th>
            <th>Number of Samples</th>
          </tr>
        </head>
        <tbody>", sep="", file=outfile)

  for (l in unique(group.labels))
  {
    cat("
          <tr>
            <td>", l, "</td>
            <td>",  length(which(group.labels == l)), "</td>
          </tr>", sep="", file=outfile)
  }

  cat("
        </tbody>
      </table>

      <h1>Sample Summary Sheets</h1>

      <p>
        For each sample a report sheet is created which summarizes the most
        relevant information using the global and local perspective.
        The global summary shows the ranked list of differentially expressed
        genes for the whole sample, the ranked list of over- and underexpressed
        gene sets after GSZ-overexpression analysis and the respective p-value
        distributions.
        The gene and gene set list are provided as tables.
      </p>

      <table>
        <thead>
          <tr>
            <th>Sample Name</th>
            <th>Group</th>
            <th>Summary Sheet</th>
            <th>Global Gene List</th>
            <th>Gene Set List</th>
          </tr>
        </thead>
        <tbody>", sep="", file=outfile)

  for (m in 1:ncol(indata))
  {
    name <- colnames(indata)[m]
    fname <- make.names(name)

    cat("
          <tr>
            <td>", name, "</td>
            <td> ", group.labels[m], "</td>
            <td>
              <a href=\"", fname, ".pdf\" target=\"_blank\">
                PDF
              </a>
            </td>
            <td>
              <a href=\"../CSV Sheets/Gene Lists - Global/", fname, ".csv\" >
                CSV
              </a>
            </td>
            <td>
              <a href=\"../CSV Sheets/Gene Set Lists - Global/", fname, ".csv\" target=\"_blank\">
                CSV
              </a>
            </td>
          </tr>", sep="", file=outfile)
  }

  cat("
        </tbody>
      </table>
    </div>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
