dri.smooth.cghdata <-
function( DNA.data, Chr, mw = 5 ) {
	DNA.smoothed <- matrix(nrow=dim(DNA.data)[1], ncol=dim(DNA.data)[2])
	for( i in c(1:24,'X','Y') ) {
		if( sum(Chr==i)>0 ) {
			m <- DNA.data[Chr==i,]
			m <- apply(m, 2, FUN=runningmean, k=mw)
			DNA.smoothed[Chr==i,] <- m
		}
	}
	return(DNA.smoothed)
}

