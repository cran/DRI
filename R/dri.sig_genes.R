dri.sig_genes <-
function(cutoff, observed, null_dist, gene_id, gene_name, chr, nuc, bt=TRUE, method="drcorrelate") {
	if(method == "drsam") {
		observed.dna <- observed$test.dna
		observed.rna <- observed$test.rna
		observed <- observed$test.summed
	
		dna.score.p <- vector()
		rna.score.p <- vector()
		summed.score.p <- vector()

		dna.score.n <- vector()
		rna.score.n <- vector()
		summed.score.n <- vector()
	}

	# positive genes
	row_number.p <- vector()
	correlation.p <- vector()
	rank.p <- vector()
	fdr_i.p <- vector()
	p.p <- vector()

	# negative genes
	row_number.n <- vector()
	correlation.n <- vector()
	rank.n <- vector()
	fdr_i.n <- vector()
	p.n <- vector()

	row_number <- c(1:length(observed))

	sorted_test <- sort(observed, decreasing=TRUE)
	sorted_test.rev <- sort(observed, decreasing=FALSE)
	sorted_abs_test <- sort(abs(observed), decreasing=TRUE)
	sorted_null <- apply(null_dist, MARGIN=2, FUN=sort, decreasing=TRUE)

	# FDR correction factor
	pi0 <- get.pi0(sorted_test, sorted_null)

	sig_count <- 0
	if( bt ) { # looking at both + and - correlations
		abs_obs.o <- order(abs(observed), decreasing = TRUE)
		row_number.s <- row_number[abs_obs.o]	
		abs_sorted_null <- abs(sorted_null)

		i <- 1
		i_pos <- 1
		i_neg <- 1
		while( sorted_abs_test[i] >= cutoff ) {
			sig_count <- sig_count + 1
			if( sorted_test[i_pos] >= abs(sorted_test.rev[i_neg]) ) { # if + correlation is more extreme than - correlation
				row_number.p <- append(row_number.p, row_number.s[i])
				correlation.p <- append(correlation.p, sorted_test[i_pos])
				rank.p <- append(rank.p, i)
				fdr_i.p <- append(fdr_i.p, fdr.fast(cutoff=sorted_test[i_pos], sig_count, abs_sorted_null, pi0))
				i_pos <- i_pos + 1
				i <- i + 1
			} else { # - correlation is more extreme than + correlation
				row_number.n <- append(row_number.n, row_number.s[i])
				correlation.n <- append(correlation.n, sorted_test.rev[i_neg])
				rank.n <- append(rank.n, i)
				fdr_i.n <- append(fdr_i.n, fdr.fast(cutoff=abs(sorted_test.rev[i_neg]), sig_count, abs_sorted_null, pi0))
				i_neg <- i_neg + 1
				i <- i + 1
			}
		}
		
		id.p <- gene_id[row_number.p]
		name.p <- gene_name[row_number.p]
		chr.p <- chr[row_number.p]
		nuc.p <- nuc[row_number.p]
		id.n <- gene_id[row_number.n]
		name.n <- gene_name[row_number.n]
		chr.n <- chr[row_number.n]
		nuc.n <- nuc[row_number.n]

		if(method == "drsam") {
			dna.score.p <- observed.dna[row_number.p]
			rna.score.p <- observed.rna[row_number.p]
			dna.score.n <- observed.dna[row_number.n]
			rna.score.n <- observed.rna[row_number.n]
			
			genes.p <- cbind(rank.p, row_number.p, id.p, name.p, chr.p, nuc.p, dna.score.p, rna.score.p, correlation.p, fdr_i.p)
			genes.n <- cbind(rank.n, row_number.n, id.n, name.n, chr.n, nuc.n, dna.score.n, rna.score.n, correlation.n, fdr_i.n)
		
	
		} else {
			genes.p <- cbind(rank.p, row_number.p, id.p, name.p, chr.p, nuc.p, correlation.p, fdr_i.p)
			genes.n <- cbind(rank.n, row_number.n, id.n, name.n, chr.n, nuc.n, correlation.n, fdr_i.n)	
		}

		Results.SigGenes <- list(positive=genes.p, negative=genes.n, cutoff=cutoff)

	} else { # only looking at + correlations
		obs.o <- order(observed, decreasing = TRUE)
		row_number.s <- row_number[obs.o]

		i <- 1
		while (observed[i] >= cutoff){
			sig_count <- sig_count + 1
			row_number.p <- append(row_number.p, row_number.s[i])
			correlation.p<- append(correlation.p, observed[i])
			rank.p <- append(rank.p, i)
			fdr_i.p <- append(fdr_i.p, fdr.fast(cutoff=observed[i], sig_count, sorted_null, pi0))
			#p <- append(p.p, pval(observed[i], null))
			i <- i + 1
		}

		id.p <- gene_id[row_number.p]
		name.p <- gene_name[row_number.p]
		chr.p <- chr[row_number.p]
		nuc.p <- nuc[row_number.p]

		genes.p <- cbind(rank.p, row_number.p, id.p, name.p, chr.p, nuc.p, correlation.p, fdr_i.p)	
		Results.SigGenes <- list(positive=genes.p, cutoff=cutoff)
	}

	return(Results.SigGenes)	
}

