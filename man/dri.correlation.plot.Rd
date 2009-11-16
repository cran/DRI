\name{dri.correlation.plot}
\alias{dri.correlation.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot significant correlations from DR-Correlate analysis}
\description{
  A plot is generated of all the correlations computed by DR-Correlate, indicating significant correlations
	in red (positive) or blue (negative), ordered by chromosomal coordinates.
}
\usage{
dri.correlation.plot(observed, Results.SigGenes, sig_cutoff, chr, 
nuc_pos, bothtails)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observed}{ vector of observed correlations from \code{drcorrelate} }
  \item{Results.SigGenes}{ list of significant genes from \code{dri.sig_genes} }
  \item{sig_cutoff}{ correlation significance cutoff returned from \code{dri.fdrCutoff} }
  \item{chr}{ vector of gene chromosome locations }
  \item{nuc_pos}{ vector of gene nucleotide positions }
  \item{bothtails}{ \code{TRUE} or \code{FALSE} indicating whether 2-tail test was performed }
}
\details{
}
\value{
  \item{plot }{ a plot of the correlations is returned in chromosomal order, with significant correlations
			marked in red (positive) and blue (negative) }
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

# Optional correlation plot for significant DR-Correlation genes
dri.correlation.plot(obs, Results, cutoff, Chr, Nuc, bothtails=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
