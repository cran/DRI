dri.merge.CNbyRNA <-
function(dna.chr, dna.nuc, dna.data, rna.chr, rna.nuc) {
	dna.nuc <- as.numeric(dna.nuc)	
	rna.nuc <- as.numeric(rna.nuc)

	dna.data <- as.matrix(dna.data)
	n <- length(dna.chr)
	average.dna.data <- matrix(nrow=length(rna.chr), ncol=dim(dna.data)[2]) # matrix full of NA
	
	i <- 1
	
	for( k in 1:length(rna.chr)) {
		queryChr <- rna.chr[k]
		queryNuc <- rna.nuc[k]

		found.indices <- FALSE
		while(!found.indices) {
			dna.chr.i <- dna.chr[i]
			dna.nuc.i <- dna.nuc[i]
			
			# RNA probe is telomeric to CGH probe on p arm
			if((dna.chr.i == queryChr) & (queryNuc < dna.nuc.i) & (i == 1 | (dna.chr[max(i-1,1)] != queryChr))) { 
				found.indices <- TRUE
				average.dna.data[k,] <- dna.data[i,] # use most proximal CGH probe available
			} else {
				
				# RNA probe is telomeric to CGH probe on q arm
				if((dna.chr.i == queryChr) & (queryNuc > dna.nuc.i) & (i == n | (dna.chr[min(i+1,n)] != queryChr))) { 
					found.indices <- TRUE
					average.dna.data[k,] <- dna.data[i,] # use most distal CGH probe available
				} else {
					
					# RNA probe is at same location as CGH probe
					if((dna.chr.i == queryChr) & (dna.nuc.i == queryNuc)) {
						found.indices <- TRUE
						average.dna.data[k,] <- dna.data[i,]
					} else {
			
						# RNA probe is between two CGH probes - most common case
						if((dna.chr.i == queryChr) & (dna.nuc.i < queryNuc) & (dna.nuc[i+1] > queryNuc)) {
							found.indices <- TRUE
							# average values from two CGH probes in each sample
							average.dna.data[k,] <- apply(dna.data[i:(i+1),], MARGIN=2, FUN=mean, na.rm=T)
						} else {
							i <- i + 1	
						}
					}
				}
			}
		}
	}
	
	average.dna.data[is.nan(average.dna.data)] <- NA # replace any NaN value with NA (occurs if both values were NA)
	
	return(average.dna.data)
}

