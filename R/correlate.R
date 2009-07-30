`correlate` <-
function(dna, rna) {
	n_genes = dim(dna)[1]
	correlation_obs <- vector(mode="numeric", length=n_genes)

	X.c <- t(scale(t(dna), center=TRUE, scale=FALSE)) # mean-centered DNA values by gene for centered Pearson correlation
	#X.c <- dna # no gene-centering for DNA data
	Y.c <- t(scale(t(rna), center=TRUE, scale=FALSE)) # mean-centered RNA values by gene for centered Pearson correlation
	#Y.c <- rna

	correlation_obs <- rowSums(X.c*Y.c, na.rm=TRUE)/(sqrt(rowSums(X.c^2, na.rm=TRUE))*sqrt(rowSums(Y.c^2, na.rm=TRUE)))

	return(correlation_obs)
}

