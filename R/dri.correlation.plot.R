`dri.correlation.plot` <-
function(observed, Results.SigGenes, sig_cutoff, chr, nuc_pos, bothtails) {
	
	# find boundaries between chrms to draw lines between them
	chr_borders <- vector() # to hold indices of first gene of each chrm
	for(i in 2:length(chr)) {
		if( chr[i] != chr[i-1] ) { # at a boundary
			chr_borders <- append(chr_borders, i)
		}
	}

	PosGenes <- Results.SigGenes$positive
	PosGenesInd <- PosGenes[,2]
	PosGenesScore <- PosGenes[,7]
	if( bothtails == TRUE ) {
		NegGenes <- Results.SigGenes$negative
		NegGenesInd <- NegGenes[,2]
		NegGenesScore <- NegGenes[,7]
	}

	par(xaxt="n", xaxs="i")
	plot(observed)
	points(x=PosGenesInd, y=PosGenesScore, col="red", pch=20)
	abline(h=sig_cutoff, lty=2, col="black")
	if( bothtails == TRUE ) {
		points(x=NegGenesInd, y=NegGenesScore, col="blue", pch=20)
		abline(h=(-1*sig_cutoff), lty=2, col="black")
	}
}

