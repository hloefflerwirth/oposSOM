pipeline.htmlSampleSummary <- function()
{
  filename <- file.path(paste(files.name, "- Results"), "Summary Sheets - Samples", "0verview.html")
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
      <h1>Error Model</h1>

      <dl>
        <dt>Error Estimation Model</dt>", sep="", file=outfile)

  if (preferences$error.model == "replicates")
  {
    cat("
        <dd>Shrinkage: SD(replicates) vs. LPE</dd>",
        sep="", file = outfile)
  } else if (preferences$error.model == "all.samples")
  {
    cat("
        <dd>LPE: SD(all genes) vs. &lt;e&gt;(all genes)</dd>",
        sep="", file = outfile)
  } else if (preferences$error.model == "groups")
  {
    cat("
        <dd>Shrinkage: SD(categories) vs. LPE</dd>",
        sep="", file = outfile)
  }


  cat("
      </dl>

      <ul>
        <li>
          <a href=\"../LPE/Sigma_LPE.pdf\" target=\"_blank\">
              LPE Error Plot For All Samples
          </a>
        </li>
      </ul>

      <h1>Categories</h1>

      <table>
        <thead>
          <tr>
            <th>Category</th>
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
        genes together with the associated significance characteristics for
        the whole sample, the ranked list of over- and underexpressed gene
        sets after GSZ-overexpression analysis and the respective p-value
        distributions.
        The local summary sheets present the analogous information for each
        single spot which is selected using the 98%-quantile criterion. The
        two maps in the left part of the sheet show the respective first
        level SOM and the selected spot, respectively.
        The full global and local lists can be downloaded in excel format.
      </p>

      <table>
        <thead>
          <tr>
            <th>Sample Name</th>
            <th>Category</th>
            <th>Summary Sheets</th>
            <th>Global Gene List</th>
            <th>Local Gene List</th>
          </tr>
        </thead>
        <tbody>", sep="", file=outfile)

  for (m in 1:ncol(indata))
  {
    cat("
          <tr>
            <td>", colnames(indata)[m], "</td>
            <td> ", group.labels[m], "</td>
            <td>
              <a href=\"", colnames(indata)[m], ".pdf\" target=\"_blank\">
                PDF
              </a>
            </td>
            <td>
              <a href=\"../CSV Sheets/Gene Lists - Global/", colnames(indata)[m], ".csv\" >
                CSV
              </a>
            </td>
            <td>", sep="", file=outfile)

    for (spot.i in 1:length(GS.infos.samples[[m]]$spots))
    {
      cat("
              <a href=\"../CSV Sheets/Gene Lists - Local/", colnames(indata)[m], ".", spot.i, ".csv\">
                CSV ", spot.i, "
              </a>", sep="", file=outfile)
    }

    cat("
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
