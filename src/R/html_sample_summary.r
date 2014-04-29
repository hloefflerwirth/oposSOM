pipeline.htmlSampleSummary <- function()
{
  filename <- file.path(paste(files.name, "- Results"), "Summary Sheets - Samples", "0verview.html")
  util.info("Writing:", filename)
  outfile <- file(filename, "w")

  cat("<!DOCTYPE html>
<html>
  <head>
    <title>Sample Summary of ", files.name, "dataset</title>
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
    <h1>Error Model</h1>

    <table>
      <tr>", sep="", file=outfile)

  if (preferences$error.model == "replicates")
  {
    cat("
        <td>Error estimation model:</td>
        <td>Shrinkage: SD(replicates) vs. LPE</td>",
        sep="", file = outfile)
  } else if (preferences$error.model == "all.samples")
  {
    cat("
        <td>Error estimation model:</td>
        <td>LPE: SD(all genes) vs. <e>(all genes)</td>",
        sep="", file = outfile)
  } else if (preferences$error.model == "groups")
  {
    cat("
        <td>Error estimation model:</td>
        <td>Shrinkage: SD(categories) vs. LPE</td>",
        sep="", file = outfile)
  }


  cat("
      </tr>
      <tr>
        <td>LPE error plot for all samples:</td>
        <td><a href=\"../LPE/Sigma_LPE.pdf\" target=\"_blank\">PDF</a></td>
      </tr>
    </table>

    <h1>Sample Summary Sheets</h1>

    <table>
      <tr>
        <td class=\"justy\">
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
        </td>
      </tr>
    </table>

    <table>
      <thead>
        <tr>
          <th style=\"width: 22%\">Sample name</th>
          <th style=\"width: 18%\">Category</th>
          <th style=\"width: 12%\">Summary Sheets</th>
          <th style=\"width: 12%\">Global Gene List</th>
          <th style=\"width: 36%\">Local Gene Lists</th>
        </tr>
      </thead>", sep="", file=outfile)

  for (m in 1:ncol(indata))
  {

    cat("
        <tr>
          <td>", colnames(indata)[m], "</td>
          <td>", group.labels[m], "</td>
          <td><a href=\"", colnames(indata)[m], ".pdf\" target=\"_blank\">PDF</a></td>
          <td><a href=\"../CSV Sheets/Gene Lists - Global/", colnames(indata)[m], ".csv\" >CSV</a></td>
          <td>", sep="", file=outfile)


    for (spot.i in 1:length(GS.infos.samples[[m]]$spots))
    {
      cat("<a href=\"../CSV Sheets/Gene Lists - Local/", colnames(indata)[m],
          ".", spot.i, ".csv\">CSV ",spot.i,"</a>", sep="", file=outfile)
    }

    cat("</td>
        </tr>", sep="", file=outfile)
  }

  cat("
      </tbody>
    </table>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
