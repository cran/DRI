`dri.fdrCutoff` <-
function(observed, null, targetFDR, bt=TRUE) {
	n_genes <- length(observed)
	
	# Take the observed/null coefficients and sort in descending order
	sorted_test <- sort(observed, decreasing=TRUE)
	sorted_null <- apply(null, MARGIN=2, FUN=sort, decreasing=TRUE)
	sorted_abs_test <- sort(abs(observed), decreasing=TRUE)

	q <- as.numeric(targetFDR)
	bt <- as.logical(bt)

	# FDR correction factor
	pi0 <- get.pi0(sorted_test, sorted_null)

	# binary search for exact number of genes with FDR < q
	if( bt ) {
		sorted_test <- sorted_abs_test
		sorted_null <- abs(sorted_null)
	} 
	n.cutoff <- binary.fdr.search(floor(n_genes/2), 0, n_genes, targetFDR=q, sorted_test, sorted_null, pi0)
	
	# return 2 element vector (n, cutoff)
	return(n.cutoff)
}

