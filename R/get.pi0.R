`get.pi0` <-
function(sorted_observed, sorted_null) {
	N <- length(sorted_observed) # number of genes in dataset
	
	# pi0 correction factor
	q75_index <- round(0.25*N)
	q25_index <- round(0.75*N)

	q75_cutoff <- median(sorted_null[q75_index,])
	q25_cutoff <- median(sorted_null[q25_index,])

	q25_75_obs <- sum(sorted_observed>q25_cutoff & sorted_observed<q75_cutoff)
	pi0 = q25_75_obs / (0.5*N)
	
	return(pi0)
}

