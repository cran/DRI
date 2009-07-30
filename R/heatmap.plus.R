`heatmap.plus` <-
function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
    distfun = dist, hclustfun = hclust, reorderfun = function(d, 
        w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
        "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
    margins = c(5, 5), ColSideColors, RowSideColorsLeft, RowSideColorsRight, cexRow = 0.2 + 
        1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, breaks, col, colsep=0, rowsep, sepcolor="white", main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
    verbose = getOption("verbose"),...) 
{
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) 
                x
            else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
        4)

    ### From heatmap2
    ## Set up breaks and force values outside the range into the endmost bins
    if(missing(breaks) || is.null(breaks) || length(breaks)<1 )
      if(missing(col))
        breaks <- 16
      else
        breaks <- length(col)+1
    if(length(breaks)==1)
      {
        breaks <- seq( min(x), max(x), length=breaks)
      }

    nbr <- length(breaks)
    ncol <- length(breaks)-1
    
    if(class(col)=="function")
      col <- col(ncol)
    else if(is.character(col) && length(col)==1)
      col <- do.call(col,list(ncol))
    
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)

    x[] <- ifelse(x<min.breaks, min.breaks, x)
    x[] <- ifelse(x>max.breaks, max.breaks, x)



    ###


    if (!missing(ColSideColors)) {
        if (!is.matrix(ColSideColors))
            stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || dim(ColSideColors)[1] != nc) 
            stop("'ColSideColors' dim()[2] must be of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColorsLeft)) {
        if (!is.matrix(RowSideColorsLeft))
            stop("'RowSideColorsLeft' must be a matrix")
        if (!is.character(RowSideColorsLeft) || dim(RowSideColorsLeft)[1] != nr) 
            stop("'RowSideColorsLeft' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    if (!missing(RowSideColorsRight)) {
        if (!is.matrix(RowSideColorsRight))
            stop("'RowSideColorsRight' must be a matrix")
        if (!is.character(RowSideColorsRight) || dim(RowSideColorsRight)[1] != nr) 
            stop("'RowSideColorsRight' must be a character vector of length nrow(x)")
        # add space for bar
        lmat <- cbind(lmat[, 1] + 1, lmat[,2] + 1, lmat[,3] + 1, c(rep(NA, nrow(lmat) - 1), 1))
	  # add space for legend
    	  lmat <- rbind(lmat[1, ] + 2, lmat[2, ] + 2, c(rep(NA, ncol(lmat)-2),1,NA), c(rep(NA, ncol(lmat)-2),2,NA))
        lwid <- append(lwid, 0.2)
    } else {
	  lmat <- rbind(lmat[1, ] + 2, lmat[2, ] + 2, c(rep(NA, ncol(lmat)-1),1), c(rep(NA, ncol(lmat)-1),2))
    }
    lhei <- append(lhei, 0.2)
    lhei <- append(lhei, 0.2)

    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, 
            "; lmat=\n")
        print(lmat)
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    # draw CGH/expr legend
    par(mar = c(2, 25, 0, 25), cex=0.75)
    min.raw <- min(x, na.rm=TRUE)
    max.raw <- max(x, na.rm=TRUE)
       
    #z <- seq(min.raw,max.raw,length=length(col))
    z <- c(-3.5,-2.5,-1.5,-.5,0,.5,1.5,2.5,3.5)
    image(z=matrix(z, ncol=1), col=col, breaks=breaks, axes=FALSE)

    par(usr=c(0,1,0,1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at=xv, labels=lv, line=-1, tick=FALSE, lty=0, cex.axis=0.75)

    # draw correlation legend
    par(mar = c(2, 25, 0, 25), cex=0.75)
          
    z <- seq(0.05, 0.95, length=10)
    corr.breaks <- seq(0,1, length=11)
    corrColorRamp.func <- colorRamp(colors=c("white","blue"), space="rgb", interpolate="spline")
    col.rgb <- corrColorRamp.func(z)
    corrCol <- vector()
    for(i in 1:10) {
	corrCol <- append(corrCol, rgb(col.rgb[i,1],col.rgb[i,2],col.rgb[i,3], max=255))
    }
    image(z=matrix(z, ncol=1), col=corrCol, breaks=corr.breaks, axes=FALSE)

    par(usr=c(0,1,0,1))
    lv <- pretty(corr.breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at=lv, labels=lv, line=-1, tick=FALSE, lty=0, cex.axis=0.75)


    # right rowsidebar goes first
    if (!missing(RowSideColorsRight)) {
        par(mar = c(1, 0.5, margins[1], 0.5))
        rsc=RowSideColorsRight;
        rsc.colors=matrix();
        rsc.names=names(table(rsc));
        rsc.i=1;
        for(rsc.name in rsc.names){
          rsc.colors[rsc.i]=rsc.name;
          rsc[rsc==rsc.name]=rsc.i;
          rsc.i=rsc.i+1;
        }
        rsc=matrix(as.numeric(rsc), nrow=dim(rsc)[1]);
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)

        if (length(colnames(RowSideColorsRight))>0) {
          	if(dim(rsc)[2] == 1) {
			axis(3, 0:(dim(rsc)[2]-1), colnames(RowSideColorsRight), las=2, line = -0.5, tick=FALSE);
		} else {
			axis(3, 0:(dim(rsc)[2]-1) / (dim(rsc)[2]-1), colnames(RowSideColorsRight), las=2, line = -0.5, tick=FALSE);
		}
        }
    }
    # left rowsidebar goes second
    if (!missing(RowSideColorsLeft)) {
        par(mar = c(1, 1.5, margins[1], 0.5))
        rsc=RowSideColorsLeft;
        rsc.colors=matrix();
        rsc.names=names(table(rsc));
        rsc.i=1;
        for(rsc.name in rsc.names){
          rsc.colors[rsc.i]=rsc.name;
          rsc[rsc==rsc.name]=rsc.i;
          rsc.i=rsc.i+1;
        }
        rsc=matrix(as.numeric(rsc), nrow=dim(rsc)[1]);
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)

        if (length(colnames(RowSideColorsLeft))>0) {
		if(dim(rsc)[2] == 1) {
			axis(3, 0:(dim(rsc)[2]-1), colnames(RowSideColorsLeft), las=2, line = -0.5, tick=FALSE);
		} else {
			axis(3, 0:(dim(rsc)[2]-1) / (dim(rsc)[2]-1), colnames(RowSideColorsLeft), las=2, line = -0.5, tick=FALSE);
		}
        }
	  if (length(rownames(RowSideColorsLeft))>0) {
		axis(2, 0:(dim(rsc)[1]-1) / (dim(rsc)[1]-1), rownames(RowSideColorsLeft), las=2, line = -0.5, tick=FALSE);
	  }
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc=ColSideColors;
        csc.colors=matrix();
        csc.names=names(table(csc));
        csc.i=1;
        for(csc.name in csc.names){
          csc.colors[csc.i]=csc.name;
          csc[csc==csc.name]=csc.i;
          csc.i=csc.i+1;
        }
        csc=matrix(as.numeric(csc), nrow=dim(csc)[1]);
        image(csc, col = as.vector(csc.colors), axes = FALSE)

        if (length(colnames(ColSideColors))>0) {
          axis(2, 0:(dim(csc)[2]-1) / (dim(csc)[2]-1), colnames(ColSideColors), las=2, tick=FALSE);
        }
    }
    #par(mar = c(margins[1], 0, margins[1], margins[2]))
    par(mar = c(1, 0, margins[1], 0))
    if (!symm || scale != "none") {
        x <- t(x)
    }
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col=col, breaks=breaks, ...)
    # axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
    axis(3, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
    
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))

    ####
    ## add 'background' colored spaces to visually separate sections
    if(!missing(colsep))
      for(csep in colsep)
        rect(xleft=csep+0.5, ybottom=rep(0,length(csep)),
             xright=csep+0.53, ytop=rep(ncol(x)+1,csep),
             lty=1, lwd=1, col=sepcolor, border=sepcolor)

    if(!missing(rowsep))
      for(rsep in rowsep)
        rect(xleft=0, ybottom=nrow(x)+1-rsep-0.5,
             xright=ncol(x)+1, ytop=nrow(x)+1-rsep-0.53,
             lty=1, lwd=1, col=sepcolor, border=sepcolor)

    ####
  
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend) 
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
    if (doCdend) 
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main)) 
        frame()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}

