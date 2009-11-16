drsam <-
function(DNA.data, RNA.data, labels, transform.type, for.null=FALSE) {
	summed.statistic <- vector()
	
	DNA.ttest <- ttest.func(x=DNA.data,y=labels)
	DNA.tt <- DNA.ttest$tt

	RNA.ttest <- ttest.func(x=RNA.data,y=labels)
	RNA.tt <- RNA.ttest$tt

	
	if( transform.type == "standardize" ) {
		z.dna <- (DNA.tt - mean(DNA.tt))/sqrt(var(DNA.tt))
		z.rna <- (RNA.tt - mean(RNA.tt))/sqrt(var(RNA.tt))		
		summed.statistic <- z.dna + z.rna
		test.dna <- z.dna # return values
		test.rna <- z.rna 
	}
	if ( transform.type == "rank" ) { # use rank
		DNA.tt.rank <- rank(DNA.tt)
		RNA.tt.rank <- rank(RNA.tt)
		summed.statistic <- (DNA.tt.rank - median(DNA.tt.rank)) + (RNA.tt.rank - median(RNA.tt.rank))
		test.dna <- (DNA.tt.rank - median(DNA.tt.rank)) # return values
		test.rna <- (RNA.tt.rank - median(RNA.tt.rank))
	}
	if( transform.type == "raw" ) { # use raw t-score
		summed.statistic <- DNA.tt + RNA.tt
		test.dna <- DNA.tt # return values
		test.rna <- RNA.tt
	}

	# weighting so genes with equal scores from DNA and RNA are favored
	ratios <- cbind(abs(test.dna/test.rna), abs(test.rna/test.dna))
	score <- apply(ratios, 1, min) * summed.statistic

	if( for.null ) {
		return(score)
	} else {
		return(list(test.summed=signif(score, digits=4), test.dna=signif(test.dna, digits=4), test.rna=signif(test.rna,digits=4)))
	}
}

