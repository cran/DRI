`permute` <-
function(dna, rna, sample_ids, k=50) {
	n_genes <- dim(dna)[1]
	n_samples <- dim(dna)[2]
	correlation_null <- matrix(nrow=n_genes, ncol=k)

	X.c <- t(scale(t(dna), center=TRUE, scale=FALSE)) # mean-centered DNA values by gene
	#X.c <- dna # no gene-centering for DNA data
	Y.c <- t(scale(t(rna), center=TRUE, scale=FALSE)) # mean-centered RNA values for permuted matrix by gene
	#Y.c <- RNA.permuted

	for(i in 1:k) {
		col_temp <- sample(sample_ids, n_samples, replace=FALSE)
		Y.c.permuted <- Y.c[,col_temp]
		
		correlation_null[,i] <- rowSums(X.c*Y.c.permuted, na.rm=TRUE)/(sqrt(rowSums(X.c^2, na.rm=TRUE))*sqrt(rowSums(Y.c.permuted^2, na.rm=TRUE)))
	}
	return(correlation_null)
}

