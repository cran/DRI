drcorrelate <-
function(DNA, RNA, method="pearson", tail_p=10) {
	DNA.data <-as.matrix(DNA)
	RNA.data <- as.matrix(RNA)

	if(method == "pearson") {
		observed <- correlate(DNA.data, RNA.data)
	}
	if(method == "spearman") {
		# convert data values to ranks then run normal pearson procedure
		DNA.data.ranks <- t(apply(DNA.data, MARGIN=1, FUN=rank))
		RNA.data.ranks <- t(apply(RNA.data, MARGIN=1, FUN=rank))

		observed <- correlate(DNA.data.ranks, RNA.data.ranks)
	}
	if(method == "ttest") {
		observed <- observed_t(DNA.data, RNA.data, tail_p)
	}	

	return(signif(observed, digits=4))
}

