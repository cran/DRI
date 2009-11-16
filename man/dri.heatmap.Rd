\name{dri.heatmap}
\alias{dri.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate a heatmap of significant DR-Correlate genes }
\description{
  A heatmap is generated showing the DNA copy number and gene expression data for the significant,
	positively correlated DR-Correlate genes.
}
\usage{
dri.heatmap(Results.SigGenes, DNA, RNA, SampleIDs, GeneNames, Chr, 
Nuc, statistic, color.scheme)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Results.SigGenes}{ list of significant genes from \code{dri.sig_genes} }
  \item{DNA}{ matrix of DNA copy number data }
  \item{RNA}{ matrix of gene expression data, samples (columns) in same order as DNA matrix }
  \item{SampleIDs}{ vector of sample names }
  \item{GeneNames}{ vector of gene names }
  \item{Chr}{ vector of gene chromosome locations }
  \item{Nuc}{ vector of gene nucleotide positions }
  \item{statistic}{ method used in \code{drcorrelate}, either "pearson", "spearman", or "ttest" }
  \item{color.scheme}{ desired heatmap color scheme, either "RG" (red-green), "RB" ("red-blue), 
				or "YB" (yellow-blue) }
}
\details{
}
\value{
}
\references{
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/
}
\author{ Keyan Salari, Robert Tibshirani, Jonathan R. Pollack }
\note{
}
\seealso{
\code{\link{drcorrelate}}, \code{\link{drcorrelate.null}}, \code{\link{drsam}}, 
\code{\link{drsam.null}}, \code{\link{dri.fdrCutoff}}, \code{\link{dri.sig_genes}}, 
\code{\link{dri.heatmap}}, \code{\link{dri.merge.CNbyRNA}}, \code{\link{dri.smooth.cghdata}}, 
\code{\link{runFusedLasso}}
}
\examples{
require(impute)
data(mySampleData)
attach(mySampleData)

# DNA data should contain no missing values - pre-smooth beforehand
# Impute missing values for gene expression data
RNA.data <- dri.impute(RNA.data)

# DR-Correlate analysis to find genes with correlated DNA/RNA measurements
obs <- drcorrelate(DNA.data, RNA.data, method="pearson")
# generate null distribution for FDR calculation (10 permutations)
null <- drcorrelate.null(DNA.data, RNA.data, method="pearson", perm=10)
# identify the correlation cutoff corresponding to your desired FDR
n.cutoff <- dri.fdrCutoff(obs, null, targetFDR=0.05, bt=TRUE)
cutoff <- n.cutoff[2]
# retrieve all genes that are significant at the determined cutoff, and
# calculate gene-specific FDRs
Results <- dri.sig_genes(cutoff, obs, null, GeneIDs, GeneNames, Chr, Nuc, 
bt=TRUE, method="drcorrelate") 

# Optional heatmap plot for significant DR-Correlation genes
sample.names <- colnames(DNA.data)
pdf(file="DRI-Heatmap.pdf", height=8, width=11)
dri.heatmap(Results, DNA.data, RNA.data, sample.names, GeneNames, Chr, Nuc, 
statistic="pearson", color.scheme="RG")
dev.off()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
