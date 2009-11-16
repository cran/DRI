fdr.fast <-
function(cutoff, n_genes, sorted_null, pi0=1) {	
	# Calculate FPs - only looks at genes more extreme than cutoff
	fp <- apply(sorted_null, MARGIN=2, FUN=sum.fp, cutoff=cutoff)
	
	FDR <- (median(fp) / n_genes) * pi0
	return(signif(FDR, digits=4))
}

