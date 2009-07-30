`permute_t` <-
function(DNA.data, RNA.data, sample_ids, tail_p, k) {
	n_genes <- dim(DNA.data)[1]
	n_samples <- dim(DNA.data)[2]
	t_null <- matrix(nrow=n_genes, ncol=k) # matrix of t-statistics for each gene and k permutations

	if( tail_p > 1 ) {
		tail_p <- tail_p/100
	}

	RNA.permuted <- matrix(nrow=n_genes, ncol=n_samples)		
	n <- ceiling(n_samples * tail_p) # this is the number of samples to take from the top and bottom of the DNA distribution

	DNA.data.o <- t(apply(DNA.data, MARGIN=1, FUN=order))
	
	for(i in 1:k) {
		col_temp <- sample(sample_ids, n_samples, replace=FALSE) # permute column labels
		
		RNA.permuted <- RNA.data[,col_temp] # construct RNA data matrix with permuted columns
		
		RNA.data.s <- RNA.data
		for( j in 1:n_genes) {
			order_indices <- DNA.data.o[j,]
			RNA.data.s[j,] <- RNA.permuted[j,order_indices]
		}

		bottom <- RNA.data.s[,1:n]
		top <- RNA.data.s[,(n_samples - n):n_samples]
	
		x <- cbind(bottom, top)
		y <- c(rep(1,n), rep(2,n))

		t_null[,i] <- ttest.func(x,y)$tt # t-statistic for t-test betw RNA values from top and bottom p of DNA data
	}

	return(t_null)
}

