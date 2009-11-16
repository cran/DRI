dri.heatmap <-
function(Results.SigGenes, DNA, RNA, SampleIDs, GeneNames, Chr, Nuc, statistic, color.scheme) {
	# statistic = c("pearson", "spearman", "ttest")
	# color.scheme = c("RG", "RB", "YB")
	
	n_genes <- dim(DNA)[1]
	n_col <- dim(DNA)[2]
	n_samples <- n_col

		#genes.p <- cbind(rank.p, row_number.p, id.p, name.p, chr.p, nuc.p, correlation.p, fdr_i.p)
		#genes.n <- cbind(rank.n, row_number.n, id.n, name.n, chr.n, nuc.n, correlation.n, fdr_i.n)	
		#Results.SigGenes <- list(positive=genes.p, negative=genes.n, cutoff=C)

	SigGenesPos <- Results.SigGenes$positive
	#RSG$positive = genes.p <- cbind(rank.p, row_number.p, id.p, name.p, correlation.p, fdr_i.p)	
	RowNumberPos <- as.numeric(SigGenesPos[,2])
	NamesPos <- SigGenesPos[,4]
	CorrelationsPos <- as.numeric(SigGenesPos[,7])

	ChrPos <- Chr[RowNumberPos]
	#ChrPos <- vector(length=length(RowNumberPos))
	#for(i in 1:length(RowNumberPos)) {
	#	ChrPos[i] <- Chr[RowNumberPos[i]]
	#}

	# positively correlated, significant genes
	physical.chr.order <- order(RowNumberPos, decreasing=FALSE)
	physical.chr.sorted <- sort(RowNumberPos, decreasing=FALSE)
	ChrPos.sorted <- ChrPos[physical.chr.order]

	# Alternate chromosome coloring with 2 colors
	ChrColorLabels <- vector(mode="character")
	last <- "plum4"
	for(i in c(1:24, "X", "Y")) {
		instances <- sum(ChrPos.sorted==i)
		if(instances > 0) {
			if(last == "plum4") {
				ChrColorLabels <- append(ChrColorLabels, rep("salmon", instances))
				last <- "salmon"
			} else {
				ChrColorLabels <- append(ChrColorLabels, rep("plum4", instances))
				last <- "plum4"
			}
		}
	}

	# chr text labels
	ChrTextLabels <- vector(mode="character")
	for(i in c(1:24, "X", "Y", "")) {
		instances <- sum(ChrPos.sorted==i)
		if(instances > 0) {
			if( instances %% 2 == 0 ) {
				midposition <- instances/2
				label <- c(rep("", midposition-1), i, rep("",midposition))
				ChrTextLabels <- append(ChrTextLabels, label)
			} else {
				midposition <- ceiling(instances/2)
				label <- c(rep("", midposition-1), i, rep("",midposition-1))
				ChrTextLabels <- append(ChrTextLabels, label)
			}
		}
	}

	corrColorRamp.WB.func <- colorRamp(colors=c("white","blue"), space="rgb", interpolate="spline")
	CorrColorLabels <- vector(mode="character", length=length(CorrelationsPos))
	
	### Decide whether to transform correlation values to 0-1 range or leave them as raw scores
	if( statistic == "ttest" ) {
		CorrelationsPos <- scale01(CorrelationsPos)
	}
	for(i in 1:length(CorrelationsPos)) {
		corr.rgb <- corrColorRamp.WB.func(CorrelationsPos[physical.chr.order[i]])
		corr.color <- rgb(corr.rgb[1],corr.rgb[2],corr.rgb[3], max=255)
		CorrColorLabels[i] <- corr.color
	}


	#index <- genes[,1]
	#dna_sample_ids <- colnames(mat_dna)[offset:n_col]
	#rna_sample_ids <- colnames(mat_rna)[offset:n_col]
	dna_sample_ids <- SampleIDs
	rna_sample_ids <- SampleIDs	

	#dna_rna_obj <- matrix(nrow=length(RowNumberPos), ncol=2*n_samples+1, dimnames=list(ChrPos[physical.chr.order], c(dna_sample_ids, "", rna_sample_ids)))
	dna_rna_obj <- matrix(nrow=length(RowNumberPos), ncol=2*n_samples, dimnames=list(ChrPos.sorted, c(dna_sample_ids, rna_sample_ids)))

	# this will order rows of dna_rna_obj in chrm order if data from Excel was in chrm order
	for(i in 1:length(RowNumberPos)) {
		dna_rna_obj[i,1:n_samples] <- as.double(DNA[physical.chr.sorted[i],])
		#dna_rna_obj[i,n_samples+1] <- 10000
		#dna_rna_obj[i,(n_samples+2):(2*n_samples+1)] <- as.double(RNA[physical.chr.sorted[i],])
		dna_rna_obj[i,(n_samples+1):(2*n_samples)] <- as.double(RNA[physical.chr.sorted[i],])
	}

	dna_rna_obj[is.na(dna_rna_obj)]<-0

	
	if( color.scheme == "RG" ) {
		col <- c(rgb(0, 4:1, 0, max=4),rgb(0,0,0,max=255),rgb(1:4, 0, 0, max=4))
	}
	if( color.scheme == "RB" ) {
		Blue2White.rgb <- corrColorRamp.WB.func(c(1,0.75,0.5,0.25))
		Blue2White <- rgb(Blue2White.rgb[,1],Blue2White.rgb[,2],Blue2White.rgb[,3], max=255)
		corrColorRamp.WR.func <- colorRamp(colors=c("white","red"), space="rgb", interpolate="spline")
		White2Red.rgb <- corrColorRamp.WR.func(c(0.25,0.5,0.75,1))
		White2Red <- rgb(White2Red.rgb[,1],White2Red.rgb[,2],White2Red.rgb[,3], max=255)
		col <- c(Blue2White,rgb(255,255,255,max=255),White2Red)
	}
	if( color.scheme == "YB" ) {
		col <- c(rgb(0, 0, 4:1, max=4),rgb(0,0,0,max=255),rgb(1:4,1:4,0, max=4))
	}

	breaks <- c(-4,-3,-2,-1,-0.25,0.25,1,2,3,4)
	rsc.left <- matrix(ChrColorLabels[length(RowNumberPos):1],nrow=length(RowNumberPos),ncol=1, dimnames=list(rev(ChrTextLabels),"Chr")) 
	rsc.right <- matrix(CorrColorLabels[length(RowNumberPos):1],nrow=length(RowNumberPos),ncol=1, dimnames=list(rep(1,length(RowNumberPos)),"Correlation"))
	heatmap.plus(dna_rna_obj[length(RowNumberPos):1,],breaks=breaks,col=col,scale='none', Rowv=NA, Colv=NA, colsep=n_samples, RowSideColorsLeft=rsc.left, RowSideColorsRight=rsc.right, labRow=NA)

}

