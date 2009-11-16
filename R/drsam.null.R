drsam.null <-
function(DNA.data, RNA.data, labels, transform.type, k) {
	n_genes <- dim(DNA.data)[1]
	n_samples <- dim(DNA.data)[2]
	tt.null <- matrix(nrow=n_genes, ncol=k) # matrix of t-statistics for each gene and k permutations

	for(i in 1:k) {
		labels.permuted <- sample(labels, n_samples, replace=FALSE) # permute sample labels
		tt.null[,i] <- drsam(DNA.data, RNA.data, labels.permuted, transform.type, for.null=TRUE)	
	}

	return(tt.null)
}

