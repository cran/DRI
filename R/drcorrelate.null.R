drcorrelate.null <-
function(DNA, RNA, method="pearson", tail_p=10, perm) {
	DNA.data <-as.matrix(DNA)
	RNA.data <- as.matrix(RNA)	
	ids <- seq(1,dim(DNA.data)[2])

	if(method == "pearson") {
		null <- permute(dna=DNA.data, rna=RNA.data, sample_ids = ids, k=as.numeric(perm))
	}
	if(method == "spearman") {
		# convert data values to ranks then run normal pearson procedure
		DNA.data.ranks <- t(apply(DNA.data, MARGIN=1, FUN=rank))
		RNA.data.ranks <- t(apply(RNA.data, MARGIN=1, FUN=rank))

		null <- permute(dna=DNA.data.ranks, rna=RNA.data.ranks, sample_ids = ids, k=as.numeric(perm))
	}
	if(method == "ttest") {
		null <- permute_t(DNA.data, RNA.data, ids, tail_p, as.numeric(perm))
	}

	return(signif(null, digits=4))
}

