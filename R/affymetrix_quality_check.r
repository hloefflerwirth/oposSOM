pipeline.affymetrixQualityCheck <- function(env)
{
  spike.in.controls = c("AFFX-LysX-3_at","AFFX-PheX-3_at","AFFX-ThrX-3_at","AFFX-DapX-3_at")
  
  if( all( spike.in.controls %in% rownames(env$indata) ) )
  {  
    
    for( i in spike.in.controls )
    {
      par(mfrow=c(2,1))
      par(mar=c(5,3,3,2))
      
      barplot( env$indata[i,], col=env$group.colors, main=paste("Intensity of sample preparation spike-in control",i), names.arg=if(ncol(env$indata)<80) colnames(env$indata) else rep("",ncol(env$indata)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(env$indata)<80,"black",NA), ylim=range( env$indata[i,] )*c(0.99,1.01), xpd=FALSE )	
      box()
      
      if( length( unique(env$group.labels) ) > 1 )
      {
        mean.boxes <- by( env$indata[i,], env$group.labels, c )[ unique( env$group.labels ) ]
        par(mar=c(5,3,0,2))
        boxplot( mean.boxes, col=env$groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n" )
        axis( 1, 1:length(env$groupwise.group.colors), unique(env$group.labels), las=2, cex.axis=0.8 )
      }
    }
    
  }
  
  
  
  if( all( c("AFFX-HUMGAPDH/M33197_3_at","AFFX-HUMGAPDH/M33197_5_at") %in% rownames(env$indata) ) )
  {  
    ratio.3.5.GAPDH <- env$indata["AFFX-HUMGAPDH/M33197_3_at",] / env$indata["AFFX-HUMGAPDH/M33197_5_at",]
    
    par(mfrow=c(2,1))
    par(mar=c(5,3,3,2))
    
    barplot( ratio.3.5.GAPDH, col=env$group.colors, main="3' / 5' ratio GAPDH", names.arg=if(ncol(env$indata)<80) colnames(env$indata) else rep("",ncol(env$indata)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(env$indata)<80,"black",NA), ylim=c( min(ratio.3.5.GAPDH )*0.99 , max(1.26,max(ratio.3.5.GAPDH )*1.01) ), xpd=FALSE )	
    abline(h=1.25, lwd=2, col="red3")
    box()
    
    if( length( unique(env$group.labels) ) > 1 )
    {
      mean.boxes <- by( ratio.3.5.GAPDH, env$group.labels, c )[ unique( env$group.labels ) ]
      par(mar=c(5,3,0,2))
      boxplot( mean.boxes, col=env$groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n" )
      axis( 1, 1:length(env$groupwise.group.colors), unique(env$group.labels), las=2, cex.axis=0.8 )
      abline(h=1.25, lwd=2, col="red3")
    }
    
  }
  
  
  
  if( all( c("AFFX-HSAC07/X00351_3_at","AFFX-HSAC07/X00351_5_at") %in% rownames(env$indata) ) )
  {	
    ratio.3.5.betaact <- env$indata["AFFX-HSAC07/X00351_3_at",] / env$indata["AFFX-HSAC07/X00351_5_at",]
    
    par(mfrow=c(2,1))
    par(mar=c(5,3,3,2))
    
    barplot( ratio.3.5.betaact, col=env$group.colors, main="3' / 5' ratio beta-actin", names.arg=if(ncol(env$indata)<80) colnames(env$indata) else rep("",ncol(env$indata)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(env$indata)<80,"black",NA), ylim=c( min(ratio.3.5.betaact )*0.99 , max(3.01,max(ratio.3.5.betaact )*1.01) ), xpd=FALSE )	
    abline(h=3, lwd=2, col="red3")
    box()
    
    if( length( unique(env$group.labels) ) > 1 )
    {
      mean.boxes <- by( ratio.3.5.betaact, env$group.labels, c )[ unique( env$group.labels ) ]
      par(mar=c(5,3,0,2))
      boxplot( mean.boxes, col=env$groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n" )
      axis( 1, 1:length(env$groupwise.group.colors), unique(env$group.labels), las=2, cex.axis=0.8 )
      abline(h=3, lwd=2, col="red3")
    }
    
  }
  
   
  
  if( all( c("AFFX-BioB-3_at","AFFX-BioC-3_at","AFFX-BioDn-3_at","AFFX-CreX-3_at") %in% rownames(env$indata) ) )
  {  
      
    par(mfrow=c(1,1))
    par(mar=c(3,3,3,3))
    plot( 0, type="n", xlim=c(1,4), ylim=range(c(env$indata["AFFX-BioB-3_at",],
                                                 env$indata["AFFX-BioC-3_at",],
                                                 env$indata["AFFX-BioDn-3_at",],
                                                 env$indata["AFFX-CreX-3_at",] )), xlab="", ylab="", axes=FALSE, main="Intensities of spike-in controls" )
    
    axis(2,las=2)
    axis(1,c(1:4),c("bioB","bioC","bioD","creX"))
    box()
    
    for( m in seq(ncol(env$indata)) )
    {
      lines( c( env$indata["AFFX-BioB-3_at",m],
                env$indata["AFFX-BioC-3_at",m],
                env$indata["AFFX-BioDn-3_at",m],
                env$indata["AFFX-CreX-3_at",m] ),
             col=env$group.colors[m] )
    }
    
  }

}

