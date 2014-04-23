
  


  ## Heatmap Wrapper cause bad Parameter-Handling Row/ColSideColors
  heatmap.wrap = function( RowSideColors, ColSideColors, ... )
  {
    if( !missing( RowSideColors ) && !missing( ColSideColors ) )
    {
      if( length( unique( RowSideColors ) ) > 1 && length( unique( ColSideColors ) ) > 1 )
      {
        heatmap( RowSideColors=RowSideColors, ColSideColors=ColSideColors, ... )
        return( NULL )
      }
    }else
    if( !missing( ColSideColors ) )
    {
      if( length( unique( ColSideColors ) ) > 1 )
      {
        heatmap( ColSideColors=ColSideColors, ... )
        return( NULL )
      }
    }  
    
    heatmap( ... )
  }




  # Filter Metagenes

  metagene.filter.list = list()


  metagene.filter.list[[1]] = list( s = c(1:preferences$dim.som1^2), n = paste("all", preferences$dim.som1^2, "Metagenes") )

  spot.metagenes = unique( unlist( sapply(GS.infos.overexpression$spots, function(x) x$metagenes ) ) )
  metagene.filter.list[[2]] = list( s = spot.metagenes, n = paste( length(spot.metagenes), "Spot-Metagenes") )
  

  if( preferences$dim.som1^2 > 1000 ) 
    metagene.filter.list[[ length(metagene.filter.list) + 1 ]] = list( s = order( apply( abs(metadata), 1, max ), decreasing=T )[1:1000], n = "High Expression: 1000 Metagenes" )
  if( preferences$dim.som1^2 > 100 ) 
    metagene.filter.list[[ length(metagene.filter.list) + 1 ]] = list( s = order( apply( abs(metadata), 1, max ), decreasing=T )[1:100], n = "High Expression: 100 Metagenes" )


  if( preferences$dim.som1^2 > 1000 ) 
    metagene.filter.list[[ length(metagene.filter.list) + 1 ]] = list( s = order( apply( abs(metadata), 1, var ), decreasing=T )[1:1000], n = "Variance: 1000 Metagenes" )
  if( preferences$dim.som1^2 > 100 ) 
    metagene.filter.list[[ length(metagene.filter.list) + 1 ]] = list( s = order( apply( abs(metadata), 1, var ), decreasing=T )[1:100], n = "Variance: 100 Metagenes" )













  pdf( paste( files.name, "- Results/2nd lvl Metagene Analysis/Similarity Analysis.pdf" ), 29.7/2.54, 21/2.54 )


  for( i in 1:length(metagene.filter.list) )
  {

    s = metadata[ metagene.filter.list[[i]]$s , ]
    par( mar=c(1,1,1,1) ) 


    
    if( ncol(metadata) > 2 )
    {  
      phylo.tree = nj( dist( t( s ) ) )
      phylo.tree$tip.label = colnames(indata)
    
      plot.phylo( phylo.tree, "unrooted", cex=0.5, tip.color=group.colors )
        title(metagene.filter.list[[i]]$n)
        box()
    }
    

    heatmap.wrap( x=s, col=colramp(1000), main=metagene.filter.list[[i]]$n, margins = c(10, 5), scale="n", labRow=NA, ColSideColors=group.colors )
      par(new=T)
      plot(0,type="n", axes=F, xlab="", ylab="" )
      legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )  

    heatmap.wrap( x=s, col=colramp(1000), main=metagene.filter.list[[i]]$n, margins = c(10, 5), scale="n", labRow=NA, ColSideColors=group.colors, Colv=NA )
      par(new=T)
      plot(0,type="n", axes=F, xlab="", ylab="" )
      legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )  


  }


  dev.off()







