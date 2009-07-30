`observed_t` <-
function(DNA.data, RNA.data, tail_p) {

	n_genes <- dim(DNA.data)[1]
	n_samples <- dim(DNA.data)[2]
	t_obs <- vector(mode="numeric", length=n_genes) # vector of t-statistics for each gene

	if( tail_p > 1 ) {
		tail_p <- tail_p/100
	}
		
	n <- ceiling(n_samples * tail_p) # this is the number of samples to take from the top and bottom of the DNA distribution
	
	DNA.data.o <- t(apply(DNA.data, MARGIN=1, FUN=order))

	RNA.data.s <- RNA.data
	for( i in 1:n_genes) {
		order_indices <- DNA.data.o[i,]
		RNA.data.s[i,] <- RNA.data[i,order_indices]
	}

	bottom <- RNA.data.s[,1:n]
	top <- RNA.data.s[,(n_samples - n + 1):n_samples]
	
	x <- cbind(bottom, top)
	y <- c(rep(1,n), rep(2,n))

	t_obs <- ttest.func(x,y)$tt

	return(t_obs)
}

