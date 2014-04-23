

smoothScatter = function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("white", 
                                                                                                    blues9)), nrpoints = 100, pch = ".", cex = 1, col = "black", 
                                   transformation = function(x) x^0.25, postPlotHook = box, 
                                   xlab = NULL, ylab = NULL, xlim, ylim, xaxs = par("xaxs"), 
                                   yaxs = par("yaxs"), ...) 
{
  if (!is.numeric(nrpoints) | (nrpoints < 0) | (length(nrpoints) != 
    1)) 
    stop("'nrpoints' should be numeric scalar with value >= 0.")
  xlabel <- if (!missing(x)) 
    deparse(substitute(x))
  ylabel <- if (!missing(y)) 
    deparse(substitute(y))
  xy <- xy.coords(x, y, xlabel, ylabel)
  xlab <- if (is.null(xlab)) 
    xy$xlab
  else xlab
  ylab <- if (is.null(ylab)) 
    xy$ylab
  else ylab
  x <- cbind(xy$x, xy$y)[is.finite(xy$x) & is.finite(xy$y), 
                         , drop = FALSE]
  if (!missing(xlim)) {
    stopifnot(is.numeric(xlim), length(xlim) == 2, is.finite(xlim))
    x <- x[min(xlim) <= x[, 1] & x[, 1] <= max(xlim), ]
  }
  else {
    xlim <- range(x[, 1])
  }
  if (!missing(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim) == 2, is.finite(ylim))
    x <- x[min(ylim) <= x[, 2] & x[, 2] <= max(ylim), ]
  }
  else {
    ylim <- range(x[, 2])
  }
  map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth, list(xlim,ylim) )
  xm <- map$x1
  ym <- map$x2
  dens <- map$fhat
  dens[] <- transformation(dens)
  dens[which(dens<0.1)]=NA    
  col = paste( colramp(256), "90", sep="" )
  image(xm, ym, z = dens, col = col, xlab = xlab, 
        ylab = ylab, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs, ...)
  if (!is.null(postPlotHook)) 
    postPlotHook()
  if (nrpoints > 0) {
    nrpoints <- min(nrow(x), ceiling(nrpoints))
    stopifnot((nx <- length(xm)) == nrow(dens), (ny <- length(ym)) == 
      ncol(dens))
    ixm <- 1L + as.integer((nx - 1) * (x[, 1] - xm[1])/(xm[nx] - 
      xm[1]))
    iym <- 1L + as.integer((ny - 1) * (x[, 2] - ym[1])/(ym[ny] - 
      ym[1]))
    sel <- order(dens[cbind(ixm, iym)])[seq_len(nrpoints)]
    points(x[sel, ], pch = pch, cex = cex, col = col)
  }
}

