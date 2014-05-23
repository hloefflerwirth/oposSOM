pipeline.2ndLvlCorrelationAnalysis <- function()
{
  bootstrap.col <- group.colors

  bootstrap.col[which(group.bootstrap.score < 80)] <-
    apply((col2rgb(bootstrap.col[which(group.bootstrap.score < 80)]) + 0.7 *
           (255 - col2rgb(bootstrap.col[which(group.bootstrap.score < 80)]))) / 255,
          2, function(x) rgb(x[1],x[2],x[3]))

  bootstrap.col[which(group.bootstrap.score < 50)] <-
    apply((col2rgb(bootstrap.col[which(group.bootstrap.score < 50)]) + 0.9 *
           (255 - col2rgb(bootstrap.col[which(group.bootstrap.score < 50)]))) / 255,
          2, function(x) rgb(x[1],x[2],x[3]))

  filename <- file.path(paste(files.name, "- Results"), "2nd lvl Metagene Analysis", "Correlation Analysis.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)

  for (i in 1:length(metagene.filter.list))
  {
    s <- metadata[metagene.filter.list[[i]]$s ,]
    par(mar=c(1, 1, 1, 1))

    # Maximum Correlation Tree
    adj.matrix <- cor(s) * -1
    g <- graph.adjacency(adj.matrix, weighted=T, mode="undirected")
    stg <- minimum.spanning.tree(g)
    layout <- layout.kamada.kawai(stg)

    plot(stg, layout=layout, vertex.size=5, vertex.label = colnames(indata),
         vertex.label.cex=if (ncol(indata)<100) 1.2 else 0.6,
         vertex.color=group.colors, main=metagene.filter.list[[i]]$n)

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")

    box()

    plot(stg, layout=layout, vertex.size=5, vertex.label = rep("",ncol(indata)),
         vertex.label.cex=if (ncol(indata)<100) 1.2 else 0.6,
         vertex.color=group.colors, main=metagene.filter.list[[i]]$n)

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")

    box()

    plot(stg, layout=layout, vertex.size=5, vertex.label = rep("",ncol(indata)),
         vertex.label.cex=if (ncol(indata)<100) 1.2 else 0.6,
         vertex.color=bootstrap.col, main=metagene.filter.list[[i]]$n)

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")

    box()

    # Collelation Backbone
    adj.matrix <- cor(s)
    diag(adj.matrix) <- 0

    adj.matrix <- apply(adj.matrix, 1, function(x)
    {
      x[order(x,decreasing=T)[-c(1:2)]] <- 0
      return(x)
    })

    adj.matrix[which(adj.matrix < 0.5)] <- 0

    g <- graph.adjacency(adj.matrix, weighted=T,  mode="undirected")
    layout <- layout.fruchterman.reingold(g)

    plot(g, layout=layout, vertex.size=5, vertex.label = colnames(indata),
         vertex.label.cex=if (ncol(indata)<100) 1.2 else 0.6,
         vertex.color=group.colors, main=metagene.filter.list[[i]]$n)

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")

    box()

    plot(g, layout=layout, vertex.size=5, vertex.label = rep("",ncol(indata)),
         vertex.label.cex=if (ncol(indata)<100) 1.2 else 0.6,
         vertex.color=group.colors, main=metagene.filter.list[[i]]$n)

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")

    box()

    plot(g, layout=layout, vertex.size=5, vertex.label = rep("",ncol(indata)),
         vertex.label.cex=if (ncol(indata)<100) 1.2 else 0.6,
         vertex.color=bootstrap.col, main=metagene.filter.list[[i]]$n)

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")

    box()

    # Correlation Network
    adj.matrix <- cor(s)
    diag(adj.matrix) <- 0
    adj.matrix[which(adj.matrix < 0.5)] <- 0

    if (max(adj.matrix) > 0)
    {
      g <- graph.adjacency(adj.matrix, weighted=T,  mode="undirected")
      layout <- layout.fruchterman.reingold(g, start=matrix(1:(2*ncol(indata)),ncol=2), niter=1000)

      plot(g, layout=layout, vertex.size=ifelse(ncol(indata) < 250, 5, 3),
           vertex.label=colnames(indata),
           vertex.label.cex=if (ncol(indata)<100) 1.2 else 0.6,
           vertex.color=group.colors, main=metagene.filter.list[[i]]$n)

      legend("bottomright", as.character(unique(group.labels)), cex=0.5,
             text.col=groupwise.group.colors, bg="white")

      box()

      plot(g, layout=layout, vertex.size=ifelse(ncol(indata)<250, 5, 3),
           vertex.label = rep("",ncol(indata)),
           vertex.label.cex=if (ncol(indata) < 100) 1.2 else 0.6,
           vertex.color=group.colors, main=metagene.filter.list[[i]]$n)

      legend("bottomright", as.character(unique(group.labels)), cex=0.5,
             text.col=groupwise.group.colors, bg="white")

      box()

      plot(g, layout=layout, vertex.size=ifelse(ncol(indata) < 250, 5, 3),
           vertex.label = rep("",ncol(indata)),
           vertex.label.cex=if (ncol(indata)<100) 1.2 else 0.6,
           vertex.color=bootstrap.col, main=metagene.filter.list[[i]]$n)

      legend("bottomright", as.character(unique(group.labels)), cex=0.5,
             text.col=groupwise.group.colors, bg="white")

      box()
    }

    # PCM
    hcl <- hclust(dist(cor(s)))

    heatmap.wrap(x=cor(s), Rowv=as.dendrogram(hcl), Colv=as.dendrogram(hcl),
                 col=colramp(1000), scale="n", main=metagene.filter.list[[i]]$n,
                 mar=c(8,8), ColSideColors=group.colors, RowSideColors=group.colors)

    par(new=T)
    plot(0,type="n", axes=F, xlab="", ylab="")

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")

    if (i <= 2)
    {
      def.par <- par(no.readonly = TRUE)

      environment(pipeline.2ndLvlModuleCorrelation) <- environment()
      pipeline.2ndLvlModuleCorrelation(s, hcl)

      par(def.par)
    }

    heatmap.wrap(x=cor(s), Rowv=NA, Colv=NA, col=colramp(1000), scale="n",
                 main=metagene.filter.list[[i]]$n, mar=c(8,8),
                 ColSideColors=group.colors, RowSideColors=group.colors)

    par(new=T)
    plot(0,type="n", axes=F, xlab="", ylab="")

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")

    o <- unlist(sapply(unique(group.labels), function(gr)
    {
      idx <- names(group.labels)[which(group.labels == gr)]

      if (length(idx) > 1)
      {
        hc <- hclust(dist(t(s[,idx])))
        return(hc$labels[hc$order])
      }
      return(idx)
    }))

    heatmap.wrap(x=cor(s)[o,o], Rowv=NA, Colv=NA, col=colramp(1000),
                 scale="n", main=metagene.filter.list[[i]]$n, mar=c(8,8),
                 ColSideColors=group.colors[o], RowSideColors=group.colors[o])

    par(new=T)
    plot(0,type="n", axes=F, xlab="", ylab="")

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")
  }

  dev.off()
}
