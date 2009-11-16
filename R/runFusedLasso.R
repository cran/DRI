runFusedLasso <-
function( DNA.data, normal.data=NA , chr, nuc, FDR ) {
	if( !is.na(normal.data) ) { 
		# filter probes from X and Y chromosomes in case normals are male-female hybes		
		X.probecount <- sum(chr==23)
		Y.probecount <- sum(chr==24)
		n <- length(chr)
		XY.filter <- c(rep(0,n - (X.probecount+Y.probecount)), rep(1,X.probecount+Y.probecount))
		Normal.FL<-cghFLasso.ref(normal.data, chromosome=chr, filter=XY.filter)
		Disease.FL<-cghFLasso(DNA.data, chromosome=chr, nucleotide.position=nuc, FL.norm=Normal.FL, FDR=FDR, filter=XY.filter)

		autosome.filter <- c(rep(1,n - (X.probecount+Y.probecount)), rep(0,X.probecount+Y.probecount))
		Disease.XY.FL<-cghFLasso(DNA.data, chromosome=chr, nucleotide.position=nuc, FDR=FDR, filter=autosome.filter)
		
		# combine autosomal data and sex chromosome data
		data.FL <- rbind(Disease.FL$Esti.CopyN, Disease.XY.FL$Esti.CopyN)
	} else {
		Disease.FL<-cghFLasso(DNA.data, chromosome=chr, nucleotide.position=nuc, FDR=FDR)
		data.FL <- Disease.FL$Esti.CopyN
	}

	return(data.FL)
}

