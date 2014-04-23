

CountExtrema = function(v)
{
  v.1 = v[2:length(v)] - v [1:(length(v)-1)]
  changes = sign(v.1[2:length(v.1)]) != sign(v.1[1:(length(v.1)-1)])
  return(sum(changes))
}

GetFirstPeak = function(v) 
{
  v.1 = v[2:length(v)] - v [1:(length(v)-1)]
  changes = sign(v.1[2:length(v.1)]) != sign(v.1[1:(length(v.1)-1)])

  return(min(which(changes)))
}



Camel.Analysis = function(data, verbose=F, folder="Density Analysis")
{

  cat("Camel Analysis:\n|")
  for (i in 1:50) cat(" ")
  cat("|\n|")
  flush.console()




  if (verbose)
  {

    dir.create(folder, showWarnings=F)

    dens.list = list()
    for (i in 1:ncol(data)) 
      dens.list[[i]] = density(data[,i]) 
    max.y <- max(unlist(lapply(dens.list,function(elem) max(elem$y))))
         max.x <- max(unlist(lapply(dens.list,function(elem) max(elem$x))))
         min.x <- min(unlist(lapply(dens.list,function(elem) min(elem$x))))
    eps = max.y/20

    pdf(paste(folder,"/Chip Densities.pdf",sep=""), 10, 10)
         plot(dens.list[[1]],xlab="probe intensity", ylab="density",col=1,ylim=c(-eps,max.y+eps),xlim=c(min.x,max.x),main="Chip Intensities",type="n")
    for (i in 1:ncol(data))
      lines(dens.list[[i]],col=1)
    dev.off()



    outfile = file(paste(folder,"/Summary.html",sep=""), "w")
            
    cat("
      <html> <head> <TITLE>Density Analysis Summary</TITLE> </head> <body bgcolor=#FFFFFF ><H1 ALIGN=CENTER >Density Analysis Summary</H1>
      <style type=text/css> p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px }
      </style> ", file = outfile)

       }


  summary = matrix(0, 9, ncol(data))
  colnames(summary) = colnames(data)
  rownames(summary) = c("delta.E.max.N", "E.max.S", "Integral.N", "Integral.S", "R", "Center.N", "Center.S", "Width.S.50.left", "Width.S.50.rigth")



  corrected.data = data
  percent.S = data


  for (i in 1:ncol(data))
  {


    den = density(data[,i], n=nrow(data))


    adj = 0.1
    while (CountExtrema(density(data[,i], adjust=adj)$y) > 3)
    {
      adj = adj+0.1  
    }


    bimod.density = density(data[,i], adjust=adj)

    N.peak.x = bimod.density$x[GetFirstPeak(bimod.density$y)]

    N.peak.y = den$y[min(which(den$x > N.peak.x))]


    p95 = N.peak.y * 0.95
        


    E95.U = min(which(den$y > p95))
    E95.O = max(which(den$y > p95))


    if (min(den$y[E95.U:E95.O]) < p95)    # 2 peaks über quantil
    {
      mi = which.min(den$y[E95.U:E95.O]) + E95.U
      E95.O = max(which(den$y[1:mi] > p95))
    }



    sum.E.p = 0
    sum.p = 0

    for (xi in E95.U:E95.O)
    {
      sum.E.p = sum.E.p + den$x[xi] * den$y[xi]
      sum.p = sum.p + den$y[xi]    
    }
    
    E.max.N = sum.E.p / sum.p

    E.max.N = den$x[which.max(den$y[E95.U:E95.O]) + E95.U]



    den$x = den$x - E.max.N

    E.max.N.old = E.max.N
    E.max.N = 0





    E.max.i = 1
    while (den$x[E.max.i] < E.max.N)
    {  
      E.max.i = E.max.i+1
    }

    den.unspec = list(x=c(), y=c())
    den.unspec$x = den$x[1:(2*E.max.i)]
    den.unspec$y = den$y[c(1:E.max.i , E.max.i:1)]

    if (length(den.unspec$x) > length(den$x))
    {  
      den.unspec$x = den.unspec$x[1:length(den$x)]
      den.unspec$y = den.unspec$y[1:length(den$x)]
    }



    den.spec = den
    for (xi in 1:length(den.unspec$x))
    {
      den.spec$y[xi] = den.spec$y[xi] - den.unspec$y[xi]
      
      if (den.spec$y[xi] < 0)
        den.spec$y[xi] = 0
    }




    den.spec$x = den.spec$x[-(1:E.max.i)]
    den.spec$y = den.spec$y[-(1:E.max.i)]


    
    if (verbose)
    {
      dir.create(paste(folder,"/Chip ",i,sep=""), showWarnings=F)

      bmp(paste(folder,"/Chip ",i,"/densities.bmp",sep=""), 800, 800)
           plot(density(data[,i]-E.max.N.old),xlab="log(Expression)", ylab="density",col="gray",main=colnames(data)[i])
      points(den.unspec$x, den.unspec$y)
      points(den.spec$x, den.spec$y, col = "blue")
      dev.off()
            }



    p95 = quantile(den.spec$y, 0.95)
        
    E95.U = min(which(den.spec$y > p95))
    E95.O = max(which(den.spec$y > p95))




    sum.E.p = 0
    sum.p = 0

    for (xi in E95.U:E95.O)
    {
      sum.E.p = sum.E.p + den.spec$x[xi] * den.spec$y[xi]
      sum.p = sum.p + den.spec$y[xi]    
    }
    
    E.max.S = sum.E.p / sum.p

    E.max.S = den.spec$x[which.max(den.spec$y[E95.U:E95.O]) + E95.U]



    I.S = 0
    for (xi in 2:length(den.spec$x))
    {
      I.S = I.S +  den.spec$y[xi] * (den.spec$x[xi] - den.spec$x[xi-1])
    }
    I.N = 0
    for (xi in 2:length(den.unspec$x))
    {
      I.N = I.N +  den.unspec$y[xi] * (den.unspec$x[xi] - den.unspec$x[xi-1])
    }

    R = I.S / (I.S + I.N) 







    sum.E.p = 0
    sum.p = 0

    for (xi in 1:length(den.spec$x))
    {
      sum.E.p = sum.E.p + den.spec$x[xi] * den.spec$y[xi]
      sum.p = sum.p + den.spec$y[xi]    
    }
    
    center.S = sum.E.p / sum.p



    sum.E.p = 0
    sum.p = 0

    for (xi in 1:length(den.unspec$x))
    {
      sum.E.p = sum.E.p + den.unspec$x[xi] * den.unspec$y[xi]
      sum.p = sum.p + den.unspec$y[xi]    
    }
    
    center.N = sum.E.p / sum.p






    h50 = max(den.spec$y) / 2
        
    E50.U = min(which(den.spec$y > h50))
    E50.O = max(which(den.spec$y > h50))


    corrected.data[,i] = corrected.data[,i] - E.max.N.old



    mix.perc = den.unspec$y[which(den.unspec$x %in% den.spec$x)] / (den.unspec$y[which(den.unspec$x %in% den.spec$x)] + den.spec$y[which(den.spec$x %in% den.unspec$x)])
    percent.S[,i] = 1 - approx(x = den.spec$x[which(den.spec$x %in% den.unspec$x)], y = mix.perc, xout = corrected.data[,i], rule=2)$y 



    summary[1, i] = E.max.N.old
    summary[2, i] = E.max.S
    summary[3, i] = I.N
    summary[4, i] = I.S
    summary[5, i] = R
    summary[6, i] = center.N
    summary[7, i] = center.S
    summary[8, i] = E.max.S - den.spec$x[E50.U]
    summary[9, i] = den.spec$x[E50.O] - E.max.S



    if (verbose)
    {



      bmp(paste(folder,"/Chip ",i,"/densities.bmp",sep=""), 800, 800)

           plot(density(data[,i]-E.max.N.old),xlab="log(Expression)", ylab="density",col="gray",main=colnames(data)[i])
      abline(v=E.max.N, lty=2, col="gray")
      abline(v=E.max.S, lty=2, col="gray")

      points(den.unspec$x, den.unspec$y)
      points(den.spec$x, den.spec$y, col = "blue")
      points(den.spec$x[E50.U], den.spec$y[E50.U], col = "red")
      points(den.spec$x[E50.O], den.spec$y[E50.O], col = "red")
      lines(c(den.spec$x[E50.U], den.spec$x[E50.O]), c(h50, h50), col = "red")

      dev.off()


      
      write.table(cbind(den$x , den$y), file=paste(folder,"/Chip ",i,"/total.txt",sep=""), sep="\t", row.names=F, col.names=F)
      write.table(cbind(den.unspec$x , den.unspec$y), file=paste(folder,"/Chip ",i,"/N.txt",sep=""), sep="\t", row.names=F, col.names=F)
      write.table(cbind(den.spec$x , den.spec$y), file=paste(folder,"/Chip ",i,"/S.txt",sep=""), sep="\t", row.names=F, col.names=F)



      cat("<h2>Chip ",i," : ",colnames(data)[i],"</h2><CENTER> <TABLE BORDER=4>
        <TR><TD><CENTER><img SRC=\"Chip ",i,"/densities.bmp\" BORDER=0 width=400 heigth=400 ></CENTER> 
        </TD></TR></TABLE></CENTER> ", sep="", file = outfile)

      cat("<CENTER><TABLE BORDER=4>
        <TR><TD> delta log Emax N </TD><TD> log Emax S </TD><TD> Integral N </TD><TD> Integral S </TD><TD> rel. Expression </TD><TD> Center N </TD><TD> Center S </TD><TD> Width S.50 left </TD><TD> Width S.50 right </TD></TR>
        <TR><TD>",E.max.N.old,"</TD><TD>",E.max.S,"</TD><TD>",I.N,"</TD><TD>",I.S,"</TD><TD>",R,"</TD><TD>",center.N,"</TD><TD>",center.S,"</TD><TD>",summary[8, i],"</TD><TD>",summary[9, i],"</TD></TR>
        </TABLE></CENTER></CENTER><br><br><br><br><br><br> ", file = outfile)
 
    }




    out.intervals = round(seq(1, ncol(data), ncol(data) * 0.02))[-1]
    out.intervals[length(out.intervals) : 49] = out.intervals[length(out.intervals)]
    cat(paste(rep("#",length(which(out.intervals == i))), collapse=""));  flush.console()


  }



  if (verbose) 
  {
    cat("<h2>All Chips</h2><CENTER> <TABLE BORDER=4> 
      <TR><TD></TD><TD> delta log Emax N </TD><TD> log Emax S </TD><TD> Integral N </TD><TD> Integral S </TD><TD> rel. Expression </TD><TD> Center N </TD><TD> Center S </TD><TD> Width S.50 left </TD><TD> Width S.50 right </TD></TR> ", file = outfile)


    for (i in 1:ncol(data))
    {
      cat("<TR><TD>",colnames(data)[i],"</TD><TD>",summary[1,i],"</TD><TD>",summary[2,i],"</TD><TD>",summary[3,i],"</TD><TD>",summary[4,i],"</TD><TD>",summary[5,i],"</TD><TD>",summary[6,i],"</TD><TD>",summary[7,i],"</TD><TD>",summary[8, i],"</TD><TD>",summary[9, i],"</TD></TR> " , file = outfile)
    }


    cat("  </TABLE></CENTER></CENTER></CENTER> </body> </html> ", file = outfile)
    close(outfile)
  }


  cat("#|\n")
  flush.console()


  ret = list("Corrected.Data"=corrected.data, "Percent.S"=percent.S, "Chip.Summary"=summary, "N.peak"=den.unspec, "S.peak"=den.spec)


  return(ret)
}




