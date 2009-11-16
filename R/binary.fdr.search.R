binary.fdr.search <-
function(n, n.left, n.right, targetFDR, sorted_test, sorted_null, pi0) {
	cutoff <- sorted_test[n]
	cutoff.plus1 <- sorted_test[n+1]
	
	n.fdr <- fdr.fast(cutoff, n, sorted_null, pi0=pi0)
	n.plus1.fdr <- fdr.fast(cutoff.plus1, (n+1), sorted_null, pi0=pi0)
	
	if( n.fdr <= targetFDR & n.plus1.fdr > targetFDR) {
		#print(paste( n, " (FDR=", n.fdr, "), ", n+1, " (FDR=", n.plus1.fdr, ")"))
		return(c(n, cutoff))
	} else {
		if( n.fdr > targetFDR ) {
			n.right <- n
			return(binary.fdr.search(n=ceiling((n.left + n)/2), n.left=n.left, n.right=n.right, targetFDR, sorted_test, sorted_null, pi0))
		} else {
			n.left <- n
			return(binary.fdr.search(n=floor((n + n.right)/2), n.left=n.left, n.right=n.right, targetFDR, sorted_test, sorted_null, pi0))
		}
	}
}

