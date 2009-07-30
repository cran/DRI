\name{drsam}
\alias{drsam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perform supervised learning analysis between 2 sample classes for DNA/RNA differences }
\description{
A test is performed to identify genes with significant differences in both DNA copy number and gene
expression between two sample groups of interest.
}
\usage{
drsam(DNA.data, RNA.data, labels, transform.type, for.null = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DNA.data}{ matrix of DNA copy number data }
  \item{RNA.data}{ matrix of gene expression data, samples (columns) in same order as DNA matrix }
  \item{labels}{ class labels of the two comparison groups, either \code{1} or \code{2} }
  \item{transform.type}{ type of transformation to apply to data, either "standardize", "rank", or "raw" }
  \item{for.null}{ used internally by \code{drsam.null}, keep as default }
}
\details{
DR-SAM (DNA/RNA-Significance Analysis of Microarrays) performs a supervised analysis to identify genes with 
statistically significant differences in both DNA copy number and gene expression between different classes 
(e.g., tumor subtype-A vs. tumor subtype-B). The goal of this analysis is to identify genetic differences 
(CNAs) that mediate gene expression differences between two groups of interest. DR-SAM implements a modified 
Student's t-test to generate for each gene two t-scores assessing differences in DNA copy number and 
differences in gene expression. A final score is computed by first summing the copy number t-score and gene 
expression t-score, and then weighting the sum by the ratio of the two t-scores. The weight is applied to 
favor genes with strong differences in both DNA copy number and gene expression between the two classes.
}
\value{
  \item{observed}{ vector of observed DR-SAM scores for each gene }
}
\references{
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/
}
\author{ Keyan Salari, Robert Tibshirani, Jonathan R. Pollack }
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
