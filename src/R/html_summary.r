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
        <dd>", ncol(indata), "</dd>
        <dt>Number of groups</dt>
        <dd>", length(unique(group.labels)), "</dd>
        <dt>Number of genes</dt>
        <dd>", nrow(indata), "</dd>
        <dt>ID type of genes</dt>
				<dd>", ifelse( preferences$database.id.type!="", preferences$database.id.type, "not defined" ), "</dd>
        <dt>Dimension of the SOM</dt>
        <dd>", preferences$dim.1stLvlSom, " x ", preferences$dim.1stLvlSom, "</dd>
        <dt>Date</dt>
        <dd>", format(Sys.time(), "%a %b %d %X %Y %Z"), "</dd>
        <dt>Analyst</dt>
        <dd>", preferences$system.info["user"], "</dd>
        <dt>oposSOM version</dt>
        <dd>", preferences$session.info$otherPkgs$oposSOM$Version, "</dd>
      </dl>

      <h1>Results</h1>
      <h2>General</h2>

      <ul>
        <li>
          <a href=\"Data Overview/Data Distribution.pdf\" target=\"_blank\">
            Data Distribution (PDF)
          </a>
        </li>
        <li>
          <a href=\"Data Overview/Chromosome Expression.pdf\" target=\"_blank\">
            Chromosome-wise gene expression heatmap (PDF)
          </a>
        </li>
      </ul>

      <h2>SOM Analysis</h2>

      <p>
        These reports comprise the SOM portraits in standard and alternative
        color scales, as well as supporting maps and profiles which provide
        supplementary information about the 1st level SOM.
      </p>", sep="", file=outfile)

  
  if( file.exists(file.path(paste(files.name, "- Results"), "Expression Portraits.pdf")) )
  {
      cat("
        <ul>
          <li>
            <a href=\"Expression Portraits.pdf\" target=\"_blank\">
              1st Level SOM Expression Portraits (PDF)
            </a>
          </li>
          <li>
            <a href=\"Expression Portraits - alternative scales.pdf\" target=\"_blank\">
              Alternative Color Scales: absolute expression, WAD, loglogFC (PDF)
            </a>
          </li>
        </ul>", sep="", file=outfile)
  }

      cat("<ul>
        <li>
          <a href=\"Supporting Maps&amp;Profiles/Supporting Maps.pdf\" target=\"_blank\">
            Supporting Maps (PDF)
          </a>
        </li>
      </ul>

      <ul>
        <li>
          <a href=\"Supporting Maps&amp;Profiles/Entropy Profiles.pdf\" target=\"_blank\">
            Entropy Profiles (PDF)
          </a>
        </li>
        <li>
          <a href=\"Supporting Maps&amp;Profiles/Topology Profiles.pdf\" target=\"_blank\">
            Topology Profiles (PDF)
          </a>
        </li>
      </ul>
      <ul>
        <li>
          <a href=\"CSV Sheets/Gene localization.csv\" target=\"_blank\">
          Localization of the genes within the SOM (CSV)
          </a>
        </li>
      </ul>", sep="", file=outfile)
      
  if( file.exists(file.path(paste(files.name, "- Results"), "Summary Sheets - Samples", "0verview.html")) )
  {
    cat("<h2>Sample Summaries</h2>
  
        <p>
          Summary page for the individual samples.
        </p>
  
        <ul>
          <li>
            <a href=\"Summary Sheets - Samples/0verview.html\" target=\"_blank\">
              Sample Reports (HTML)
            </a>
          </li>
        </ul>", sep="", file=outfile)
  }
  
  if( file.exists(file.path(paste(files.name, "- Results"), "Geneset Analysis", "0verview.html")) )
  {
    cat("
      <h2>Geneset Enrichment Analysis</h2>

      <p>
        Functional analyses using predefined gene sets. The results are
        visualized in terms of heatmaps, profile plots and population maps.
      </p>

      <ul>
        <li>
          <a href=\"Geneset Analysis/0verview.html\" target=\"_blank\">
            Functional Analysis (HTML)
          </a>
        </li>
      </ul>", sep="", file=outfile)
  }
  
      
  if( file.exists(file.path(paste(files.name, "- Results"), "Sample Similarity Analysis")) ) 
  {      
    cat("
        <h2>Sample Similarity Analyses</h2>
  
        <p>
          Sample similarity analyses based on expression module data and metagene data using different metrics and algorithms.<br>
          <br>
          Euclidean Distance base approaches:
        </p>
  
        <ul>
          <li>
            <a href=\"Sample Similarity Analysis/Sample SOM.pdf\" target=\"_blank\">
              Sample SOM (PDF)
            </a>
          </li>
          <li>
            <a href=\"Sample Similarity Analysis/t-SNE.pdf\" target=\"_blank\">
              t-SNE dimension reduction (PDF)
            </a>
          </li>
          <li>
            <a href=\"Sample Similarity Analysis/Neighbor Joining.pdf\" target=\"_blank\">
              Phylogenetic clustering using neighbor joining algorithm (PDF)
            </a>
          </li>
          <li>
            <a href=\"Sample Similarity Analysis/Hierarchical Clustering.pdf\" target=\"_blank\">
              Expression heatmaps and hierarchical clustering dendrograms (PDF)
            </a>
          </li>
        </ul>
  
        Correlation base approaches:
  
        <ul>
          <li>
            <a href=\"Sample Similarity Analysis/Correlation Spanning Tree.pdf\" target=\"_blank\">
              Correlation spanning tree connecting all samples (PDF)
            </a>
          </li>
          <li>
            <a href=\"Sample Similarity Analysis/Correlation Backbone.pdf\" target=\"_blank\">
              Correlation backbone, a 2-nearest-neighbor graph (PDF)
            </a>
          </li>
          <li>
            <a href=\"Sample Similarity Analysis/Correlation Network.pdf\" target=\"_blank\">
              Correlation network connecting samples with r > 0.5 (PDF)
            </a>
          </li>
          <li>
            <a href=\"Sample Similarity Analysis/Correlation Maps.pdf\" target=\"_blank\">
              Pairwise correlation maps (supervised & clustered orderings) (PDF)
            </a>
          </li>
        </ul>
  
        Component Analysis:
    
        <ul>
          <li>
            <a href=\"Sample Similarity Analysis/Component Analysis.pdf\" target=\"_blank\">
              Independent Component Analysis: 2d-ICA, 3d-ICA (PDF)
            </a>
          </li>
        </ul>", sep="", file=outfile)
  }  
    
  cat("
    <h2>Expression Module Analysis</h2>

    <p>
      Different criteria of spot module definition such as overexpression or
      mutual correlations between the metagenes where applied.
      The reports comprise integrated portraits, functional analyses.
    </p>

    <ul>
      <li>
        <a href=\"Summary Sheets - Modules/0verview.html\" target=\"_blank\">
          Spot Reports (HTML)
        </a>
      </li>
    </ul>

    <p>
      Analyses based on PAT-wise aggregated data, including clustering dendrogram, 
      portraits and assotiation to the sample groups.
    </p>
      
      <ul>
        <li>
        <a href=\"Summary Sheets - PATs/PAT report.pdf\" target=\"_blank\">
        PAT Report (PDF)
        </a>
      </li>
    </ul>", sep="", file=outfile)    


  if( file.exists(file.path(paste(files.name, "- Results"), "Summary Sheets - Groups/0verview.html")) ) 
  {  
    cat("
      <h2>Group Analyses</h2>
      <p>
        Analyses based on group-wise aggregated data, including portraits,
        clustering and functional analyses.
      </p>

      <ul>
        <li>
          <a href=\"Summary Sheets - Groups/0verview.html\" target=\"_blank\">
            Group Analysis Reports (HTML)
          </a>
        </li>
      </ul>", sep="", file=outfile)
  }

  if( file.exists(file.path(paste(files.name, "- Results"), "Summary Sheets - Differences/0verview.html")) ) 
  {
    cat("
      <h2>Pairwise Differences Analyses</h2>

      <p>
        Analyses based on pairwise comparisons, including portraits,
        clustering and functional analyses.
      </p>

      <ul>
        <li>
          <a href=\"Summary Sheets - Differences/0verview.html\" target=\"_blank\">
            Differences Analysis Reports (HTML)
          </a>
        </li>
      </ul>", sep="", file=outfile)
  }

  cat("
    </div>
  </body>
</html>", sep="", file=outfile)

  close(outfile)
}
