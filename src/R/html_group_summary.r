pipeline.htmlGroupSummary <- function()
{
  if (length(unique(group.labels)) <= 1)
  {
    return()
  }

  filename <- file.path(paste(files.name, "- Results"),
                        "Summary Sheets - Groups",
                        "0verview.html")

  util.info("Writing:", filename)
  f <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Group Summary of ", files.name, " dataset</title>
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
      <h1>Group Centered Analyses</h1>

      <p>
        The SOM portraits are aggregated in groupwise fashions. Further,
        group stability scores are calculated using correlation silhouette method.
        Finally, samples within each group are clustered and visualized as
        hierarchical dendrograms.
      </p>

      <ul>
        <li>
          <a href=\"Expression Portraits Groups.pdf\" target=\"_blank\">
            Expression Portraits (PDF)
          </a>
        </li>
        <li>
          <a href=\"Group Assignment.pdf\" target=\"_blank\">
            Assignment Stability  Scores (PDF)
          </a>
        </li>
        <li>
          <a href=\"Group Clustering.pdf\" target=\"_blank\">
            Groupwise Clustering Dendrograms (PDF)
          </a>
        </li>
      </ul>

      <h2>Group Summary Sheets</h2>

      <p>
        For each group a report sheet, gene lists and gene set lists are
        created analogous to the sample summary sheets. Additionally, specific
        gene sets are given which show highest significance in statistical
        testing of the group compared to all other groups.
      </p>

      <table>
        <thead>
          <tr>
            <th>Group</th>
            <th>Summary Sheet</th>
            <th>Global Gene List</th>
            <th>Gene Set List</th>", sep="", file=f)

  if (preferences$activated.modules$geneset.analysis)
  {
    cat("
            <th colspan=\"2\">Specific Gene Set</th>", sep="", file=f)
  }

  cat("
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
            </td>", sep="", file=f)

    if (preferences$activated.modules$geneset.analysis)
    {
      cat("
            <td>
              <a href=\"Geneset Analysis/Specific GS ", fname, ".pdf\" target=\"_blank\">
                PDF
              </a>
            </td>
            <td>
              <a href=\"Geneset Analysis/Specific GS ", fname, ".csv\" target=\"_blank\">
                CSV
              </a>
            </td>", sep="", file=f)
    }

    cat("
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
