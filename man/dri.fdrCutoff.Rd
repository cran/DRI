\name{dri.fdrCutoff}
\alias{dri.fdrCutoff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determine cutoff score for a desired false discovery rate (FDR) }
\description{
  A search is performed for the cutoff score (for either DR-Correlate or DR-SAM)
	that corresponds to the user-defined FDR.
}
\usage{
dri.fdrCutoff(observed, null, targetFDR, bt = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observed}{ vector of scores from either \code{drcorrelate} or \code{drsam} }
  \item{null}{ matrix of null data from either \code{drcorrelate.null} or \code{drsam.null} }
  \item{targetFDR}{ desired false discovery rate }
  \item{bt}{ either \code{TRUE} or \code{FALSE} indicating whether a 2-tail test was performed }
}
\details{
A binary search is implemented to find the cutoff score that corresponds to the user-defined FDR
}
\value{
  \item{comp1 }{ a two element list containing the number of genes found significant at the chosen
			FDR, and the score cutoff corresponding to that FDR. }
}
\references{
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/
}
\author{ Keyan Salari, Robert Tibshirani, and Jonathan R. Pollack }
\note{
}
\seealso{
\code{\link{drcorrelate}},\code{\link{drcorrelate.null}},\code{\link{drsam}},
\code{\link{drsam.null}},\code{\link{dri.fdrCutoff}},\code{\link{dri.sig_genes}},\code{\link{dri.heatmap}}
}
\examples{
require(impute)
data(mySampleData)
attach(mySampleData)

# DNA data should contain no missing values - pre-mooth beforehand
# Impute missing values for gene expression data
RNA.data <- dri.impute(RNA.data)

# DR-Correlate analysis to find genes with correlated DNA/RNA measurements
obs <- drcorrelate(DNA.data, RNA.data, method="pearson")
# generate null distribution for FDR calculation (10 permutations)
null <- drcorrelate.null(DNA.data, RNA.data, method="pearson", perm=10)
# identify the correlation cutoff corresponding to your desired FDR
n.cutoff <- dri.fdrCutoff(obs, null, targetFDR=0.05, bt=TRUE)
cutoff <- n.cutoff[2]
# retrieve all genes that are significant at the determined cutoff, and calculate gene-specific FDRs
Results <- dri.sig_genes(cutoff, obs, null, GeneIDs, GeneNames, Chr, Nuc, bt=TRUE, method="drcorrelate") 

# DR-SAM analysis to find genes with alterations in both DNA and RNA between different classes
labels <- c(rep(1,25), rep(2,25)) # first 25 samples belong to class 1 and next 25 to class 2
obs <- drsam(DNA.data, RNA.data, labels, transform.type="raw")
# generate null distribution for FDR calculation (10 permutations)
null <- drsam.null(DNA.data, RNA.data, labels, transform.type="raw", 10)
# identify the correlation cutoff corresponding to your desired FDR
n.cutoff <- dri.fdrCutoff(obs$test.summed, null, targetFDR=0.05, bt=TRUE)
cutoff <- n.cutoff[2]
# retrieve all genes that are significant at the determined cutoff, and calculate gene-specific FDRs
Results <- dri.sig_genes(cutoff, obs, null, GeneIDs, GeneNames, Chr, Nuc, bt=TRUE, method="drsam") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
